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

#include "FemPhysicsLinearElasticity.h"
#include <iostream>
#include <Engine/Profiler.h>

#pragma warning( disable : 4267) // for size_t-uint conversions
#pragma warning( disable : 4244) // for double-float conversions

//#define USE_TRACE_FORMULA

namespace FEM_SYSTEM
{
	FemPhysicsLinearElasticity::FemPhysicsLinearElasticity(std::vector<Tet>& tetrahedra,
		std::vector<Node>& nodes, const FemConfig& config)
		: FemPhysicsLinear(tetrahedra, nodes, config)
	{
		if (mSimType == ST_EXPLICIT)
		{
			// compute decomposition for mass matrix
			mMassMatrix.makeCompressed();
			mSparseLU.compute(mMassMatrix);
			bool failed = mSparseLU.info() != Eigen::Success;
			if (failed)
			{
				Printf("LU decomposition failed\n");
				EigenMatrix denseM(mMassMatrix);
				{
					real det = denseM.determinant();
					auto fullPivLU = denseM.fullPivLu();
					auto rank = fullPivLU.rank();
					Printf("rank: %d, det: %f\n", rank, det);
				}
			}
		}
		
		AssembleStiffnessMatrix();
	}

	void FemPhysicsLinearElasticity::Step(real dt)
	{
		if (mSimType == ST_STATIC)
		{
			if (mForceFraction == 0)
			{
				mForceFraction = 1;
				SolveEquilibrium(1);
				CheckForInversion();
			}
		}
		else if (mSimType == ST_QUASI_STATIC)
		{
			if (mForceFraction < 1)
			{
				mForceFraction = std::min(real(1), mForceFraction + mForceStep);
				SolveEquilibrium(mForceFraction);
				CheckForInversion();
			}
		}
		else
		{
			mForceFraction = std::min(real(1), mForceFraction + mForceStep);
			real h = dt / mNumSteps;
			for (int i = 0; i < mNumSteps; i++)
			{
				mPreviousPositions = mDeformedPositions;
				if (mSimType == ST_EXPLICIT)
					StepExplicit(h);
				else if (mSimType == ST_IMPLICIT)
					StepImplicit(h);
				HandleCollisions(h);
			}
			CheckForInversion();
		}
	}

	void FemPhysicsLinearElasticity::SolveEquilibrium(float t)
	{
		EigenMatrix K(mStiffnessMatrix);
		auto decomp = K.lu();

		// add the gravity forces
		EigenVector f = ComputeLoadingForces();

		// add term due to Dirichlet BCs
#ifdef USE_COMPUTE_FORCES
		std::vector<Vector3R> fe(GetNumFreeNodes());
		ComputeForces(fe);
		f += GetEigenVector(fe);
#else
		f -= mForceFraction * mBCStiffnessMatrix * GetEigenVector(mInitialDisplacements);
#endif

		EigenVector sol = decomp.solve(t * f);
		Vector3Array u = GetStdVector(sol);
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			mDeformedPositions[i] = mReferencePositions[i + mNumBCs] + u[i];
		}
	}

	// stepper optimized for solving linear elasticity dynamics problems: Mu'' + Ku = f
	void FemPhysicsLinearElasticity::StepExplicit(real h)
	{
		// TODO: avoid allocations by using class data members

		// compute nodal displacements and aggregate them into one big vector
		size_t numNodes = GetNumFreeNodes();
		Vector3Array u(numNodes);
		for (size_t i = 0; i < numNodes; i++)
		{
			u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
		}

		// compute the elastic forces using the linear elasticity relation f = Ku
		// (when using explicit time integration we don't really need to assemble the K matrix - see ComputeForces(),
		// but it's good to have it for other purposes, e.g. damping, implicit integration, quasi-static simulation)
		EigenVector f = -mStiffnessMatrix * GetEigenVector(u);

		// add the gravity forces
		f += ComputeLoadingForces();

		// solve the linear system for the accelerations (invert the mass matrix) 
		EigenVector sol = mSparseLU.solve(f);
		Vector3Array a = GetStdVector(sol);

		// apply the forces and integrate (Symplectic Euler)
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] += h * a[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsLinearElasticity::StepImplicit(real h)
	{
		PROFILE_SCOPE("LE implicit");
		// compute nodal displacements and aggregate them into one big vector
		size_t numNodes = GetNumFreeNodes();
		Vector3Array u(numNodes);
		for (size_t i = 0; i < numNodes; i++)
		{
			u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
		}

		// compute the elastic forces using the linear elasticity relation f = Ku
		// (when using explicit time integration we don't really need to assemble the K matrix - see ComputeForces(),
		// but it's good to have it for other purposes, e.g. damping, implicit integration, quasi-static simulation)
		EigenVector f = -mStiffnessMatrix * GetEigenVector(u);

		// add the gravity forces
		f += ComputeLoadingForces();

		// build the RHS
		EigenVector b = mMassMatrix * GetEigenVector(mVelocities) + h * f;

		// solve the linear system
		if (!mHaveSysMatrix)
		{
			mHaveSysMatrix = true;
			if (mUseImplicitPressureForces)
				mSolver.Init(mMassMatrix + h * h * (mStiffnessMatrix + mTractionStiffnessMatrix), LST_LU, true);
			else
				mSolver.Init(mMassMatrix + h * h * mStiffnessMatrix, LST_LU, true);
		}
		EigenVector v = mSolver.Solve(b);
		Vector3Array vels = GetStdVector(v);

		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] = vels[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	// assemble the global stiffness matrix for linear elasticity (using linear elements)
	void FemPhysicsLinearElasticity::AssembleStiffnessMatrix()
	{
		size_t numNodes = GetNumFreeNodes(); // remove the first 4 entries which are fixed (hack)
		size_t numDofs = numNodes * 3;
		mStiffnessMatrix.resize(numDofs, numDofs);
		mStiffnessMatrix.setZero();
		mBCStiffnessMatrix.resize(numDofs, mNumBCs * 3);
		mBCStiffnessMatrix.setZero();
		// go through all linear elements (tetrahedra)
		for (size_t i = 0; i < (size_t)mTetMesh->GetNumElements(); i++)
		{
			EigenMatrix Klocal;
			ComputeLocalStiffnessMatrixBB(i, GetShearModulus(), GetLameFirstParam(), Klocal);
			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (size_t j = 0; j < GetNumLocalNodes(); j++)
			{
				size_t jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, j)];
				for (size_t k = 0; k < GetNumLocalNodes(); k++)
				{
					size_t kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, k)];
					if (jGlobal < mNumBCs)
						continue;
					int jOffset = (jGlobal - mNumBCs) * NUM_POS_COMPONENTS;
					int kOffset = (kGlobal - mNumBCs) * NUM_POS_COMPONENTS;

					// add the the whole 3x3 block to the block matrix
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							real val = Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
							if (kGlobal < mNumBCs)
							{
								mBCStiffnessMatrix.coeffRef(jOffset + x, kGlobal * NUM_POS_COMPONENTS + y) += val;
							}
							else
							{
								mStiffnessMatrix.coeffRef(jOffset + x, kOffset + y) += val;
							}
						}
					}
				}
			}
		}
	}

	void FemPhysicsLinearElasticity::ComputeLocalStiffnessMatrix(uint32 i, EigenMatrix& Klocal)
	{
		// this code only works for linear elements
		ASSERT(GetNumLocalNodes() == 4);
		Klocal = EigenMatrix(12, 12);

		Matrix3R Bn[4], Bs[4];
		ComputeStrainJacobian(i, Bn, Bs);

		// for linear FEM we can actually precompute the tangent stiffness matrix
		real s = mYoungsModulus / (1.f + mPoissonRatio);
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				Matrix3R block = mElementVolumes[i] * (Bn[j] * mNormalElasticityMatrix * Bn[k] + s * Bs[j] * Bs[k]);
				// copy the block into Klocal
				for (uint32 a = 0; a < NUM_POS_COMPONENTS; a++)
				{
					for (uint32 b = 0; b < NUM_POS_COMPONENTS; b++)
					{
						Klocal(j * NUM_POS_COMPONENTS + a, k * NUM_POS_COMPONENTS + b) = block(a, b);
					}
				}
			}
		}
	}

	void FemPhysicsLinearElasticity::ComputeLocalStiffnessMatrixBB1(uint32 k, EigenMatrix& Klocal)
	{
		// this code only works for linear elements
		ASSERT(GetNumLocalNodes() == 4);
		Klocal = EigenMatrix(12, 12);

		real s = mYoungsModulus / (1.f + mPoissonRatio);
		// convert elastic params to Lame coefficients [Sifakis]
		const real mu = 0.5f * s;
		const real lambda = s * mPoissonRatio / (1.f - 2 * mPoissonRatio);

		Matrix3R Bn[4], Bs[4];
		ComputeStrainJacobian(k, Bn, Bs); // TODO: call only once

		Matrix3R identity;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
#ifdef USE_TRACE_FORMULA
				Matrix3R L = Matrix3R::TensorProduct(mBarycentricJacobians[k].y[i], mBarycentricJacobians[k].y[j]);
				Matrix3R block = mElementVolumes[k] * (lambda * L + mu * !L + mu * L.Trace() * identity);
#else
				//Matrix3R Kn = Bn[i] * Bn[j];
				//Matrix3R Ks = Bs[i] * Bs[j];
				//Matrix3R block = mElementVolumes[k] * ((lambda + 2 * mu) * Kn + mu * Ks); // wrong!!!
				Matrix3R block = mElementVolumes[k] * (Bn[i] * mNormalElasticityMatrix * Bn[j] + mu * Bs[i] * Bs[j]);
#endif

				// copy the block into Klocal
				for (uint32 a = 0; a < NUM_POS_COMPONENTS; a++)
				{
					for (uint32 b = 0; b < NUM_POS_COMPONENTS; b++)
					{
						Klocal(i * NUM_POS_COMPONENTS + a, j * NUM_POS_COMPONENTS + b) = block(a, b);
					}
				}
			}
		}
	}

	void FemPhysicsLinearElasticity::ComputeDyadicMatrixBB2(uint32 i, uint32 j, Vector3R y[4], Matrix3R& L)
	{
		int* multiIndexI = mTetMesh->GetIJKL(i);
		int* multiIndexJ = mTetMesh->GetIJKL(j);
		// compute every component of L
		for (uint32 a = 0; a < NUM_POS_COMPONENTS; a++)
		{
			for (uint32 b = 0; b < NUM_POS_COMPONENTS; b++)
			{
				real l = 0;
				for (uint32 c = 0; c < NUM_BARYCENTRIC_COMPONENTS; c++)
				{
					for (uint32 d = 0; d < NUM_BARYCENTRIC_COMPONENTS; d++)
					{
						MultiIndex multiIndex1(multiIndexI);
						MultiIndex multiIndex2(multiIndexJ);
						if (multiIndex1[c] == 0 || multiIndex2[d] == 0)
							continue;
						multiIndex1[c]--;
						multiIndex2[d]--;
						l += y[c][a] * y[d][b] * ComputeMultiIndexSumFactor(2, multiIndex1, multiIndex2);
					}
				}
				L(a, b) = l;
			}
		}
	}

	void FemPhysicsLinearElasticity::ComputeLocalStiffnessMatrixBB2(uint32 k, EigenMatrix& Klocal)
	{
		ASSERT(GetNumLocalNodes() == 10);
		Klocal = EigenMatrix(30, 30);

		const real mu = GetShearModulus();
		const real lambda = GetLameFirstParam();

		Matrix3R Bn[4], Bs[4];
		ComputeStrainJacobian(k, Bn, Bs);

		Matrix3R identity;
		for (uint32 i = 0; i < GetNumLocalNodes(); i++)
		{
			for (uint32 j = 0; j < GetNumLocalNodes(); j++)
			{
#ifdef USE_TRACE_FORMULA
				Matrix3R L;
				ComputeDyadicMatrixBB2(i, j, mBarycentricJacobians[k].y, L);				
				Matrix3R block = 0.4f * mElementVolumes[k] * (lambda * L + mu * !L + mu * L.Trace() * identity);
#else
				const auto& y = mBarycentricJacobians[k].y;
				int* multiIndexI = mTetMesh->GetIJKL(i);
				int* multiIndexJ = mTetMesh->GetIJKL(j);
				Matrix3R Kn(0), Ks(0);
				for (uint32 c = 0; c < NUM_BARYCENTRIC_COMPONENTS; c++)
				{				
					for (uint32 d = 0; d < NUM_BARYCENTRIC_COMPONENTS; d++)
					{
						MultiIndex multiIndex1(multiIndexI);
						MultiIndex multiIndex2(multiIndexJ);
						if (multiIndex1[c] == 0 || multiIndex2[d] == 0)
						{
							continue;
						}
						multiIndex1[c]--;
						multiIndex2[d]--;
						real factor = ComputeMultiIndexSumFactor(2, multiIndex1, multiIndex2);
						Kn = Kn + factor * Bn[c] * mNormalElasticityMatrix * Bn[d];
						Ks = Ks + factor * mu * Bs[c] * Bs[d];
					}
				}
				Matrix3R block = 0.4 * mElementVolumes[k] * (Kn + Ks);
#endif

				// copy the block into Klocal
				for (uint32 a = 0; a < NUM_POS_COMPONENTS; a++)
				{
					for (uint32 b = 0; b < NUM_POS_COMPONENTS; b++)
					{
						// the factor of 2 was put in there to obtain the proper nodal forces - its origin is still unclear
						Klocal(i * NUM_POS_COMPONENTS + a, j * NUM_POS_COMPONENTS + b) = 4 * block(a, b);
					}
				}
			}
		}
	}

} // namespace FEM_SYSTEM