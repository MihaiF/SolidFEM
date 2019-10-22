/*
BSD 3-Clause License

Copyright (c) 2019, Mihai Francu
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "FemPhysicsLinear.h"
#include "MeshFactory.h"
#include <Engine/Profiler.h>
#include <iostream>

#pragma warning( disable : 4244) // for double-float conversions
#pragma warning( disable : 4267) // for size_t-uint conversions

// References:
// [Weber] Weber, D. et al., Interactive deformable models with quadratic bases in Bernstein-Bezier form

// TODO: add Rayleigh damping

namespace FEM_SYSTEM
{
	FemPhysicsLinear::FemPhysicsLinear(std::vector<Tetrahedron>& tetrahedra,
		std::vector<Node>& nodes, const FemConfig& config)
		: FemPhysicsBase(config)
		, mOrder(config.mOrder)
		, mUseLumpedMass(true)
	{
		ASSERT(tetrahedra.size() != 0);
		ASSERT(nodes.size() != 0);
		CreateMeshAndDofs(tetrahedra, nodes);

		mElementVolumes.resize(GetNumElements());
		mBarycentricJacobians.resize(GetNumElements());
		for (uint32 e = 0; e < GetNumElements(); e++)
		{
			ComputeBarycentricJacobian(e, mBarycentricJacobians[e].y);
		}

		if (mSimType > ST_QUASI_STATIC)
			AssembleMassMatrix();

		mBodyForces.resize(GetNumFreeNodes());
		ComputeBodyForces(mBodyForces);
	}

	void FemPhysicsLinear::CreateMeshAndDofs(std::vector<Tetrahedron>& tetrahedra, std::vector<Node>& nodes)
	{
#ifdef LOG_TIMES
		Printf(("----- FemPhysicsAnyOrderCorotationalElasticity for order " + std::to_string(order) + "\n").c_str());
#endif // LOG_TIMES
		// 1. Build the higher order mesh/convert to a TetrahedralMesh instance
		// use the linear mesh stored as an array of Tetrahedron in the base class
		StridedVector<uint16> stridedVec(&(tetrahedra[0].i[0]), (uint)tetrahedra.size(), sizeof(Tetrahedron));

		double mesh_time;
		{
			MEASURE_TIME_P("mesh_time", mesh_time);
			mTetMesh.reset(MeshFactory::MakeMesh<TetrahedralMesh<uint32>>(mOrder, nodes.size(), stridedVec, tetrahedra.size()));
		}

		std::vector<Vector3R> points(nodes.size());
		double other_time;
		{
			MEASURE_TIME_P("other_time", other_time);
			// 2. Compute the interpolated positions
			// copy over the node positions of the linear mesh to a linear array
			// TODO: could use a strided vector instead
			for (uint32 i = 0; i < nodes.size(); i++)
			{
				points[i] = nodes[i].pos;
			}
		}

		std::vector<Vector3R> interpolatedPositions(GetNumNodes());
		double interpolate_time;
		{
			MEASURE_TIME_P("interpolate_time", interpolate_time);
			// interpolate the mesh node positions (for the given order)
			// the first resulting nodes are the same as in the linear mesh
			MeshFactory::Interpolate(*mTetMesh, &points[0], &interpolatedPositions[0]);
		}

		std::vector<bool> fixed(GetNumNodes(), false);
		double other_time2;
		{
			MEASURE_TIME_P("other_time2", other_time2);
			// 3. Mark all boundary nodes
			// initialize the fixed flags array
			for (uint32 i = 0; i < nodes.size(); i++)
			{
				fixed[i] = nodes[i].invMass == 0;
			}
		}
		other_time += other_time2;

		double boundary_time;
		{
			MEASURE_TIME_P("boundary_time", boundary_time);
			// 3.5 Check the higher order nodes and ensure that the expected leftmost cantilever nodes are marked fixed
			if (mOrder > 1)
			{
				for (uint32 e = 0; e < GetNumElements(); e++)
				{
					// TODO These are hacky solutions, need supported iterators from the mesh class..

					// First check all edge nodes
					std::tuple<int, int> edge_pairs[]{
						std::make_tuple(0, 1),
						std::make_tuple(0, 2),
						std::make_tuple(0, 3),
						std::make_tuple(1, 2),
						std::make_tuple(1, 3),
						std::make_tuple(2, 3)
					};

					int no_nodes_pr_edge = mTetMesh->GetNodesPerEdge(mOrder);
					int edge_startidx = 4;
					for (int i = 0; i < 6; i++)
					{
						int c1, c2; std::tie(c1, c2) = edge_pairs[i];
						uint32 gidx_c1 = mTetMesh->GetGlobalIndex(e, c1);
						uint32 gidx_c2 = mTetMesh->GetGlobalIndex(e, c2);
						//if (nodes[gidx_c1].invMass == 0 && nodes[gidx_c2].invMass == 0)
						if (fixed[gidx_c1] && fixed[gidx_c2])
						{
							for (int j = 0; j < no_nodes_pr_edge; j++)
							{
								int edge_lidx = edge_startidx + i * no_nodes_pr_edge + j;
								uint32 gidx_edgej = mTetMesh->GetGlobalIndex(e, edge_lidx);
								fixed[gidx_edgej] = true;
							}
						}
					}

					// Then check all face nodes, if any
					std::tuple<int, int, int> face_pairs[]{
						std::make_tuple(0, 1, 2), std::make_tuple(0, 1, 3),
						std::make_tuple(0, 2, 3), std::make_tuple(1, 2, 3)
					};

					int no_nodes_pr_face = mTetMesh->GetNodesPerFace(mOrder);
					int face_startidx = 4 + no_nodes_pr_edge * 6;
					for (int i = 0; i < 4; i++)
					{
						int c1, c2, c3; std::tie(c1, c2, c3) = face_pairs[i];

						uint32 gidx_c1 = mTetMesh->GetGlobalIndex(e, c1);
						uint32 gidx_c2 = mTetMesh->GetGlobalIndex(e, c2);
						uint32 gidx_c3 = mTetMesh->GetGlobalIndex(e, c3);

						if (fixed[gidx_c1] && fixed[gidx_c2] && fixed[gidx_c3])
						{
							for (int j = 0; j < no_nodes_pr_face; j++)
							{
								int face_lidx = face_startidx + i * no_nodes_pr_face + j;
								uint32 gidx_facej = mTetMesh->GetGlobalIndex(e, face_lidx);
								fixed[gidx_facej] = true;
							}
						}
					}
				}
			}
		}

		// 4. Create a index-mapping from mesh-global-index to index into the mReferencePosition and mDeformedPositions vectors
		// - this is too allow all of the fixed nodes to be listed first in the two vectors.

		double other_time3;
		{
			MEASURE_TIME_P("other_time3", other_time3);
			// create a mapping with the fixed nodes first
			mReshuffleMapInv.clear();
			for (uint32 i = 0; i < fixed.size(); i++)
			{
				if (fixed[i])
					mReshuffleMapInv.push_back(i);
			}
			mNumBCs = mReshuffleMapInv.size();
			for (uint32 i = 0; i < fixed.size(); i++)
			{
				if (!fixed[i])
					mReshuffleMapInv.push_back(i);
			}
			// create the reference positions - first ones are the fixed ones
			mReferencePositions.resize(GetNumNodes());
			mReshuffleMap.resize(GetNumNodes());
			for (uint32 i = 0; i < GetNumNodes(); i++)
			{
				uint32 idx = mReshuffleMapInv[i];
				mReferencePositions[i] = interpolatedPositions[idx];
				mReshuffleMap[idx] = i;
			}

			// 5. Copy the mReferencePositions into mDeformedPositions, to initialize it
			// - note that mDeformedPositions does not contain the boundary/fixed nodes.
			// create a vector of deformed positions (dofs)
			// the first mNumBCs nodes are fixed so we only consider the remaining ones
			mDeformedPositions.resize(GetNumFreeNodes());
			mDeformedPositions.assign(mReferencePositions.begin() + mNumBCs, mReferencePositions.end());

			// 6. Initialize the velocities vector, note that it does not consider the fixed nodes either
			// velocities should be zero by construction
			mVelocities.resize(GetNumFreeNodes());
			for (uint32 i = mNumBCs; i < nodes.size(); i++)
			{
				mVelocities[i - mNumBCs] = nodes[i].vel;
			}
		}
		other_time += other_time3;
	}

	void FemPhysicsLinear::ComputeBarycentricJacobian(uint32 i, Vector3R y[4])
	{
		uint32 i0 = mTetMesh->GetGlobalIndex(i, 0);
		uint32 i1 = mTetMesh->GetGlobalIndex(i, 1);
		uint32 i2 = mTetMesh->GetGlobalIndex(i, 2);
		uint32 i3 = mTetMesh->GetGlobalIndex(i, 3);
		const Vector3R& x0 = mReferencePositions[mReshuffleMap[i0]];
		const Vector3R& x1 = mReferencePositions[mReshuffleMap[i1]];
		const Vector3R& x2 = mReferencePositions[mReshuffleMap[i2]];
		const Vector3R& x3 = mReferencePositions[mReshuffleMap[i3]];
		Vector3R d1 = x1 - x0;
		Vector3R d2 = x2 - x0;
		Vector3R d3 = x3 - x0;
		Matrix3R mat(d1, d2, d3); // this is the reference shape matrix Dm [Sifakis][Teran]
		Matrix3R X = mat.GetInverse(); // Dm^-1
		real vol = abs(mat.Determinant()) / 6.f; // volume of the tet
		mElementVolumes[i] = vol;

		// compute the gradient of the shape functions [Erleben][Mueller]
		y[1] = X[0];
		y[2] = X[1];
		y[3] = X[2];
		y[0] = y[1] + y[2] + y[3];
		y[0].Flip();
	}

	real FemPhysicsLinear::GetTotalVolume() const
	{
		real totalVol = 0;
		for (uint32 e = 0; e < GetNumElements(); e++)
		{
			uint32 i0 = mTetMesh->GetGlobalIndex(e, 0);
			uint32 i1 = mTetMesh->GetGlobalIndex(e, 1);
			uint32 i2 = mTetMesh->GetGlobalIndex(e, 2);
			uint32 i3 = mTetMesh->GetGlobalIndex(e, 3);
			const Vector3R& x0 = GetDeformedPosition(i0);
			const Vector3R& x1 = GetDeformedPosition(i1);
			const Vector3R& x2 = GetDeformedPosition(i2);
			const Vector3R& x3 = GetDeformedPosition(i3);
			Vector3R d1 = x1 - x0;
			Vector3R d2 = x2 - x0;
			Vector3R d3 = x3 - x0;
			Matrix3R mat(d1, d2, d3); // this is the spatial shape matrix Ds [Sifakis][Teran]
			real vol = abs(mat.Determinant()) / 6.f; // volume of the tet			
			totalVol += vol;
		}
		return totalVol;
	}

	void FemPhysicsLinear::ComputeStrainJacobian(uint32 i, Matrix3R Bn[4], Matrix3R Bs[4])
	{
		Vector3R* y = mBarycentricJacobians[i].y;

		// compute the Cauchy strain Jacobian matrix according to [Mueller]: B = de/dx
		for (int j = 0; j < 4; j++)
		{
			Bn[j] = Matrix3R(y[j]); // diagonal matrix
			Bs[j] = Symm(y[j]);
		}
	}

	// compute local mass matrix (using linear shape functions)
	void ComputeLocalMassMatrix(real density, uint32 numLocalNodes, EigenMatrix& Mlocal)
	{
		// this is actually wrong for (i,i) terms - see ComputeLocalMassMatrixBB1
		ASSERT(numLocalNodes == 4);
		real massDiv = density / 20.f;
		for (size_t i = 0; i < numLocalNodes; i++)
		{
			size_t offsetI = i * 3;
			for (size_t j = 0; j < numLocalNodes; j++)
			{
				size_t offsetJ = j * 3;
				// build a diagonal matrix
				for (size_t k = 0; k < 3; k++)
					Mlocal(offsetI + k, offsetJ + k) = (i == j) ? 2 * massDiv : massDiv;
			}
		}
		//std::cout << "Mlocal" << std::endl << Mlocal << std::endl;
	}

	void ComputeLocalMassMatrixLumped(real density, uint32 numLocalNodes, EigenMatrix& Mlocal)
	{
		// this is actually wrong for (i,i) terms - see ComputeLocalMassMatrixBB1
		ASSERT(numLocalNodes == 4);
		// 4x4 identity blocks times the total mass / 20
		real massDiv = density / 4.0;
		for (size_t i = 0; i < numLocalNodes * 3; i++)
		{
			Mlocal(i, i) = massDiv;			
		}
	}

	// // compute local mass matrix using any order Bernstein-Bezier shape functions
	void FemPhysicsLinear::ComputeLocalMassMatrixBB(real density, uint32 numLocalNodes, EigenMatrix& Mlocal)
	{
		int *i_ijkl, *j_mnop;
		size_t i, j, k;
		real factor = density * (1.f / Binom(2 * mOrder + 3, 3));
		for (i = 0; i < numLocalNodes; i++)
		{
			i_ijkl = mTetMesh->GetIJKL(i);
			for (j = 0; j < numLocalNodes; j++)
			{
				j_mnop = mTetMesh->GetIJKL(j);
				real mass = factor * G(i_ijkl, j_mnop, mOrder);

				for (k = 0; k < 3; k++)
				{
					Mlocal(i * 3 + k, j * 3 + k) = mass;
				}
			}
		}
	}

	// compute local mass matrix using quadratic Bernstein-Bezier shape functions
	void FemPhysicsLinear::ComputeLocalMassMatrixBB2(real density, uint32 numLocalNodes, EigenMatrix& Mlocal)
	{
		ASSERT(numLocalNodes == 10);
		for (size_t i = 0; i < numLocalNodes; i++)
		{
			size_t offsetI = i * 3;
			auto multiIndexI = mTetMesh->mIJKL + i * 4;
			for (size_t j = 0; j < numLocalNodes; j++)
			{
				//if (i != j) continue;
				auto multiIndexJ = mTetMesh->mIJKL + j * 4;
				uint32 c1 = Combinations(multiIndexI[0] + multiIndexJ[0], multiIndexI[0]);
				uint32 c2 = Combinations(multiIndexI[1] + multiIndexJ[1], multiIndexI[1]);
				uint32 c3 = Combinations(multiIndexI[2] + multiIndexJ[2], multiIndexI[2]);
				uint32 c4 = Combinations(multiIndexI[3] + multiIndexJ[3], multiIndexI[3]);
				real massDiv = c1 * c2 * c3 * c4 * density / 210.f;
				//real massDiv = density * ComputeMultiIndexSumFactor(2, multiIndexI, multiIndexJ) / 10.f;
				// build a diagonal matrix
				size_t offsetJ = j * 3;
				for (size_t k = 0; k < 3; k++)
					Mlocal(offsetI + k, offsetJ + k) = massDiv;
			}
		}
	}

	// assemble the global mass matrix 
	void FemPhysicsLinear::AssembleMassMatrix()
	{
		// compute the local mass matrix first
		size_t numLocalNodes = GetNumLocalNodes();
		size_t numLocalDofs = NUM_POS_COMPONENTS * numLocalNodes;

		// !NOTE assumed to be constant/equal for all tetrahedrons (except for the volume) -> true for linear rest shapes and subparametric formulation
		EigenMatrix Mlocal(numLocalDofs, numLocalDofs);
		Mlocal.setZero();
		if (mOrder == 1 && mUseLumpedMass)
			ComputeLocalMassMatrixLumped(mDensity, numLocalNodes, Mlocal);
		else
			ComputeLocalMassMatrixBB(mDensity, numLocalNodes, Mlocal);

		size_t numNodes = GetNumFreeNodes();
		size_t numDofs = numNodes * 3;
		EigenMatrix M(numDofs, numDofs);
		M.setZero();
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			real vol = mElementVolumes[i];
			// for every local node there corresponds a 3x3 block matrix
			// we go through all of these blocks so that we know what the current pair is
			for (size_t j = 0; j < numLocalNodes; j++)
			{
				for (size_t k = 0; k < numLocalNodes; k++)
				{
					size_t jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, j)];
					size_t kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, k)];
					// do not add to the global matrix if at least one of the nodes is fixed
					if (jGlobal < mNumBCs || kGlobal < mNumBCs)
						continue;
					// add the local 3x3 matrix to the global matrix
					int jOffset = (jGlobal - mNumBCs) * NUM_POS_COMPONENTS;
					int kOffset = (kGlobal - mNumBCs) * NUM_POS_COMPONENTS;
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							M.coeffRef(jOffset + x, kOffset + y) += vol * Mlocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
						}
					}
				}
			}
		}

		mMassMatrix = M.sparseView();
	}

	// compute the contribution of gravity (as a force distribution) to the current forces (using BB shape functions)
	void FemPhysicsLinear::ComputeBodyForces(Vector3Array& f)
	{
		Vector3R flocal = mDensity * mGravity * (1.f / Binom(mOrder + 3, 3));
		for (size_t e = 0; e < GetNumElements(); e++)
		{
			// We note that the body force is constant across the element.
			real vol = mElementVolumes[e];
			Vector3R bforce = vol * flocal;
			for (size_t i = 0; i < GetNumLocalNodes(); i++)
			{
				size_t globalI = mReshuffleMap[mTetMesh->GetGlobalIndex(e, i)];
				if (globalI < mNumBCs)
					continue;
				int offset = (globalI - mNumBCs);
				f[offset] += bforce;
			}
		}
	}

	void FemPhysicsLinear::ComputeDeformationGradient(uint32 e, Matrix3R& F)
	{
		// compute deformed/spatial shape matrix Ds [Sifakis]
		uint32 i0 = mTetMesh->GetGlobalIndex(e, 0);
		uint32 i1 = mTetMesh->GetGlobalIndex(e, 1);
		uint32 i2 = mTetMesh->GetGlobalIndex(e, 2);
		uint32 i3 = mTetMesh->GetGlobalIndex(e, 3);
		// TODO: avoid calling virtual function
		const Vector3R& x0 = GetDeformedPosition(i0);
		const Vector3R& x1 = GetDeformedPosition(i1);
		const Vector3R& x2 = GetDeformedPosition(i2);
		const Vector3R& x3 = GetDeformedPosition(i3);
		Vector3R d1 = x1 - x0;
		Vector3R d2 = x2 - x0;
		Vector3R d3 = x3 - x0;
		Matrix3R Ds(d1, d2, d3);
		Matrix3R X(mBarycentricJacobians[e].y[1], mBarycentricJacobians[e].y[2], mBarycentricJacobians[e].y[3]);
		// compute deformation gradient
		F = Ds * !X;
	}

	// computes the local stiffness matrix using Bernstein-Bezier polynomials given the Lame parameters
	void FemPhysicsLinear::ComputeLocalStiffnessMatrixBB(uint32 elementidx, real mu, real lambda, EigenMatrix& Klocal)
	{
		uint32 numLocalNodes = GetNumLocalNodes();
		uint32 NUM_POS_COMPONENTS = 3;
		Klocal = EigenMatrix(numLocalNodes * NUM_POS_COMPONENTS, numLocalNodes * NUM_POS_COMPONENTS);

		Vector3R* barycentric_over_x = mBarycentricJacobians[elementidx].y;
		real common_factor = (3.f * mOrder * mElementVolumes[elementidx]) / (4.f * mOrder * mOrder - 1.f);

		for (uint32 i = 0; i < GetNumLocalNodes(); i++)
		{
			int *i_ijkl = mTetMesh->GetIJKL(i);
			int ijkl_c[4]{ 0 }; memcpy(ijkl_c, i_ijkl, 4 * sizeof(int));
			for (uint32 j = 0; j < GetNumLocalNodes(); j++)
			{
				int *j_mnop = mTetMesh->GetIJKL(j);
				int mnop_d[4]{ 0 }; memcpy(mnop_d, j_mnop, 4 * sizeof(int));
				for (uint32 a = 0; a < NUM_POS_COMPONENTS; a++)
				{
					for (uint32 b = 0; b < NUM_POS_COMPONENTS; b++)
					{
						real sum1 = 0;
						for (int c = 0; c < 4; c++)
						{
							ijkl_c[c] -= 1;
							for (int d = 0; d < 4; d++)
							{
								mnop_d[d] -= 1;
								sum1 += barycentric_over_x[c][a] * barycentric_over_x[d][b] * ComputeMultiIndexSumFactor(mOrder - 1, ijkl_c, mnop_d);
								mnop_d[d] += 1;
							}
							ijkl_c[c] += 1;
						}

						real sum2 = 0;
						for (int c = 0; c < 4; c++)
						{
							ijkl_c[c] -= 1;
							for (int d = 0; d < 4; d++)
							{
								mnop_d[d] -= 1;
								sum2 += barycentric_over_x[c][b] * barycentric_over_x[d][a] * ComputeMultiIndexSumFactor(mOrder - 1, ijkl_c, mnop_d);
								mnop_d[d] += 1;
							}
							ijkl_c[c] += 1;
						}
						real sum3 = 0;
						if (a == b)
						{
							for (int c = 0; c < 4; c++)
							{
								ijkl_c[c] -= 1;
								for (int d = 0; d < 4; d++)
								{
									mnop_d[d] -= 1;
									real Gv = ComputeMultiIndexSumFactor(mOrder - 1, ijkl_c, mnop_d);
									for (int k = 0; k < 3; k++)
									{
										sum3 += barycentric_over_x[c][k] * barycentric_over_x[d][k] * Gv;
									}
									mnop_d[d] += 1;
								}
								ijkl_c[c] += 1;
							}
						}

						real V_ij_ab = lambda * common_factor * sum1;
						real U_ij_ab = mu * common_factor * sum2;
						real W_ij_ab = mu * common_factor * sum3;
						real K_ij_ab = V_ij_ab + U_ij_ab + W_ij_ab;

						Klocal(i * NUM_POS_COMPONENTS + a, j * NUM_POS_COMPONENTS + b) = K_ij_ab;
					}
				}
			}
		}
	}

	void FemPhysicsLinear::SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, const std::vector<uint32>& elemList, real pressure)
	{
		mAppliedPressure = pressure;
		mTractionSurface.resize(triangleList.size());

		for (uint32 i = 0; i < triangleList.size(); i++)
		{
			int origIdx = triangleList[i]; // index in the original linear mesh
			int idx = mReshuffleMap[origIdx] - mNumBCs; // index in the reshuffled displacement nodes array
			ASSERT(idx >= 0);
			mTractionSurface[i] = idx;
		}
	}

	void FemPhysicsLinear::ComputeTractionForces()
	{
		if (mTractionSurface.empty())
			return;

		if (mUseImplicitPressureForces)
		{
			uint32 numDofs = GetNumFreeNodes() * NUM_POS_COMPONENTS;
			mTractionStiffnessMatrix.resize(numDofs, numDofs);
			mTractionStiffnessMatrix.setZero();
		}
			
		// accumulate forces
		mTractionForces.resize(GetNumFreeNodes());
		for (uint32 i = 0; i < mTractionForces.size(); i++)
			mTractionForces[i].SetZero();

		// go through all inner boundary triangles
		for (uint32 t = 0; t < mTractionSurface.size() / 3; t++)
		{
			int base = t * 3;
			// global (shuffled) indices of nodes
			uint32 i1 = mTractionSurface[base];
			uint32 i2 = mTractionSurface[base + 1];
			uint32 i3 = mTractionSurface[base + 2];

			// compute triangle area and normal
			const Vector3R& p1 = mDeformedPositions[i1];
			const Vector3R& p2 = mDeformedPositions[i2];
			const Vector3R& p3 = mDeformedPositions[i3];

			auto computeForce = [&](const Vector3R& p1, const Vector3R& p2, const Vector3R& p3)->Vector3R
			{
				Vector3R a = p2 - p1;
				Vector3R b = p3 - p1;
				Vector3R normal = cross(a, b);
				real area = 0.5f * normal.Length();
				normal.Normalize();

				Vector3R traction = mAppliedPressure * normal;
				Vector3R force = (area / 3.0f) * traction;
				return force;
			};

			// apply traction to triangle nodes
			Vector3R force = computeForce(p1, p2, p3);
			mTractionForces[i1] += force;
			mTractionForces[i2] += force;
			mTractionForces[i3] += force;

			if (mUseImplicitPressureForces)
			{
				// compute local stiffness matrix
				Matrix3R K[3];
				Vector3R a = p2 - p1;
				Vector3R b = p3 - p1;
				a.Scale(mAppliedPressure / 6.0f);
				b.Scale(mAppliedPressure / 6.0f);
				K[2] = Matrix3R::Skew(a);
				K[1] = Matrix3R::Skew(-b);
				K[0] = -K[1] - K[2];

				// finite difference approximation - for verification purposes
				//real dx = 1e-4f;
				//real invDx = 1.f / dx;
				//Vector3R dfdx = invDx * (computeForce(p1, p2, p3 + Vector3R(dx, 0, 0)) - force);
				//Vector3R dfdy = invDx * (computeForce(p1, p2, p3 + Vector3R(0, dx, 0)) - force);
				//Vector3R dfdz = invDx * (computeForce(p1, p2, p3 + Vector3R(0, 0, dx)) - force);
				//Matrix3 K2(dfdx, dfdy, dfdz);

				//dfdx = invDx * (computeForce(p1, p2 + Vector3R(dx, 0, 0), p3) - force);
				//dfdy = invDx * (computeForce(p1, p2 + Vector3R(0, dx, 0), p3) - force);
				//dfdz = invDx * (computeForce(p1, p2 + Vector3R(0, 0, dx), p3) - force);
				//Matrix3 K1(dfdx, dfdy, dfdz);

				//dfdx = invDx * (computeForce(p1 + Vector3R(dx, 0, 0), p2, p3) - force);
				//dfdy = invDx * (computeForce(p1 + Vector3R(0, dx, 0), p2, p3) - force);
				//dfdz = invDx * (computeForce(p1 + Vector3R(0, 0, dx), p2, p3) - force);
				//Matrix3 K0(dfdx, dfdy, dfdz);

				// assemble global stiffness matrix
				for (uint32 x = 0; x < 3; x++)
				{
					for (uint32 y = 0; y < 3; y++)
					{
						mTractionStiffnessMatrix(i1 * 3 + x, i1 * 3 + y) += K[0](x, y);
						mTractionStiffnessMatrix(i1 * 3 + x, i2 * 3 + y) += K[1](x, y);
						mTractionStiffnessMatrix(i1 * 3 + x, i3 * 3 + y) += K[2](x, y);
						mTractionStiffnessMatrix(i2 * 3 + x, i1 * 3 + y) += K[0](x, y);
						mTractionStiffnessMatrix(i2 * 3 + x, i2 * 3 + y) += K[1](x, y);
						mTractionStiffnessMatrix(i2 * 3 + x, i3 * 3 + y) += K[2](x, y);
						mTractionStiffnessMatrix(i3 * 3 + x, i1 * 3 + y) += K[0](x, y);
						mTractionStiffnessMatrix(i3 * 3 + x, i2 * 3 + y) += K[1](x, y);
						mTractionStiffnessMatrix(i3 * 3 + x, i3 * 3 + y) += K[2](x, y);
					}
				}
			}
		}
	}

	EigenVector FemPhysicsLinear::ComputeLoadingForces()
	{
		// add the gravity forces
		EigenVector f = GetEigenVector(mBodyForces);

		// add boundary traction forces
		ComputeTractionForces();
		if (mTractionForces.size() > 0)
		{
			f += GetEigenVector(mTractionForces);
		}

		f *= mForceFraction; // slow application of forces

		return f;
	}

	bool FemPhysicsLinear::CheckForInversion()
	{
		const bool verbose = false;
		bool ret = false;
		if (verbose)
			Printf("Check for inversion\n");
		for (uint32 i = 0; i < GetNumElements(); i++)
		{
			Matrix3R F;
			ComputeDeformationGradient(i, F);
			real det = F.Determinant();
			if (verbose)
				Printf("%f\n", det);
			if (det <= 0.f)
			{
				ret = true;
				//Printf("Inverted element %d\n", i);
			}
		}

		if (ret)
			Printf("Inversion detected\n");

		return ret;
	}

} // namespace FEM_SYSTEM
