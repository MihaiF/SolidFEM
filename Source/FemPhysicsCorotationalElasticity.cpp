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

// Custom references
#include "FemPhysicsCorotationalElasticity.h"
#include "MeshFactory.h"
#include "PolarDecomposition.h"
#include <Engine/Profiler.h>

#pragma warning( disable : 4267) // for size_t-uint conversions
#pragma warning( disable : 4244) // for double-float conversions

//#define LOG_TIMES

namespace FEM_SYSTEM
{
	FemPhysicsCorotationalElasticity::FemPhysicsCorotationalElasticity(std::vector<Tet>& tetrahedra,
		std::vector<Node>& nodes, const FemConfig& config)
		: FemPhysicsLinear(tetrahedra, nodes, config)
	{
		double stiffnessmat_time;
		{
			MEASURE_TIME_P("stiffnessmat_time", stiffnessmat_time);
			CacheLocalStiffnessMatrices();
		}

		if (mSimType == ST_EXPLICIT)
			mExplicitMassSolver.Init(mMassMatrix, LST_LU);
	}

	void FemPhysicsCorotationalElasticity::Step(real dt)
	{
		// static analysis would be the same as for linear elasticity (unless we use a nonlinear corotational formulation)
		if (mSimType == ST_STATIC || mSimType == ST_QUASI_STATIC)
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
			mForceFraction = 1; // std::min(1.f, mForceFraction + mForceStep);
			real h = dt / mNumSteps;
			for (int i = 0; i < mNumSteps; i++)
			{
				mPreviousPositions = mDeformedPositions; // for collisions
				if (mSimType == ST_EXPLICIT)
					StepExplicit(h);
				else if (mSimType == ST_IMPLICIT)
					StepImplicit(h);
				HandleCollisions(h);
			}
			CheckForInversion();
		}
	}

	// stepper optimized for solving linear elasticity dynamics problems: Mu'' + Ku = f
	void FemPhysicsCorotationalElasticity::StepExplicit(real h)
	{
		PROFILE_SCOPE("CR explicit");
		size_t numNodes = GetNumFreeNodes();

		ComputeRotations();
		ComputeElasticForces(mElasticForce);

		// add the gravity forces
		EigenVector f = ComputeLoadingForces() - mElasticForce;

		// solve the linear system for the accelerations (invert the mass matrix) 
		EigenVector sol = mExplicitMassSolver.Solve(f);
		Vector3Array a = GetStdVector(sol);

		// apply the forces and integrate (Symplectic Euler)
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] += h * a[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsCorotationalElasticity::StepImplicit(real h)
	{
		PROFILE_SCOPE("AnyOrder CR implicit");
		ComputeRotations();
		AssembleStiffnessMatrix();
		auto S = mMassMatrix + h * h * mStiffnessMatrix;

		size_t numNodes = GetNumFreeNodes();

		ComputeElasticForces(mElasticForce);
		EigenVector f = ComputeLoadingForces() - mElasticForce;

		// build the RHS
		EigenVector b = mMassMatrix * GetEigenVector(mVelocities) + h * f;

		// solve the linear system
		EigenVector v;
		{
			PROFILE_SCOPE("CR solve");
#if !defined(_DEBUG) && defined(USE_MKL)
			Eigen::PardisoLU<SparseMatrix> decomp;
#else
			Eigen::SparseLU<SparseMatrix> decomp;
#endif
			decomp.compute(S);
			v = decomp.solve(b);
		}
		Vector3Array vels = GetStdVector(v);

		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] = vels[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsCorotationalElasticity::SolveEquilibrium(float t)
	{
#ifdef LOG_TIMES
		Printf("----- FemPhysicsAnyOrderCorotationalElasticity SolveEquilibrium (Quasi-static) \n");
#endif // LOG_TIMES

		double stiffnessmat_time = 0;
		double decomp_time = 0;
		double elasticf_time = 0;
		double solve_time = 0;
		double update_time = 0;

		double tmp_time;

		{
			Eigen::SimplicialLLT<SparseMatrix> decomp;

			{
				MEASURE_TIME_P("stiffnessmat_time", tmp_time);
				ComputeRotations();
				AssembleStiffnessMatrix();
			}
			stiffnessmat_time += tmp_time;

			{
				MEASURE_TIME_P("decomp_time", tmp_time);
				decomp.compute(mStiffnessMatrix);
			}
			decomp_time += tmp_time;

			EigenVector elasticForce;
			EigenVector f;
			{
				MEASURE_TIME_P("elasticf_time", tmp_time);
				ComputeElasticForces(elasticForce);
				mForceFraction = t;
				EigenVector bodyForces = ComputeLoadingForces();
				f = bodyForces - elasticForce;
			}
			elasticf_time += tmp_time;

			EigenVector sol;

			{
				MEASURE_TIME_P("solve_time", tmp_time);
				sol = decomp.solve(f);
			}
			solve_time += tmp_time;

			{
				MEASURE_TIME_P("update_time", tmp_time);
				Vector3Array u = GetStdVector(sol);
				for (uint32 i = 0; i < GetNumFreeNodes(); i++)
				{
					mDeformedPositions[i] += u[i];
				}
			}
			update_time += tmp_time;
		}

#ifdef LOG_TIMES
		double totaltime = stiffnessmat_time + decomp_time + elasticf_time + solve_time + update_time;
		Printf(("- Time for assembling stiffness matrix = " + std::to_string(stiffnessmat_time*1000.f) + "us ->" + std::to_string(stiffnessmat_time / totaltime) + "%%\n").c_str());
		Printf(("- Time for decomposition of K = " + std::to_string(decomp_time*1000.f) + "us = "+std::to_string(decomp_time/1000.f)+"s ->" + std::to_string(decomp_time / totaltime) + "%%\n").c_str());
		Printf(("- Time for computing elastif forces = " + std::to_string(elasticf_time*1000.f) + "us ->" + std::to_string(elasticf_time / totaltime) + "%%\n").c_str());
		Printf(("- Time for solving the system = " + std::to_string(solve_time*1000.f) + "us ->" + std::to_string(solve_time / totaltime) + "%%\n").c_str());
		Printf(("- Time for updating the positions = " + std::to_string(update_time*1000.f) + "us ->" + std::to_string(update_time / totaltime) + "%%\n").c_str());
		Printf(("- Totaltime = " + std::to_string(totaltime*1000.f) + "us = " + std::to_string(totaltime) + "ms = " + std::to_string(totaltime / 1000.f) + "s\n").c_str());
#endif // LOG_TIMES
	}

	void FemPhysicsCorotationalElasticity::ComputeElasticForces(EigenVector& elasticForce) const
	{
		PROFILE_SCOPE("CR compute f");
		size_t numNodes = GetNumFreeNodes(); // remove the first 4 entries which are fixed (hack)
		size_t numDofs = numNodes * 3;
		elasticForce.resize(numDofs, 1);
		elasticForce.setZero();

		// go through all elements (tetrahedra)
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			EigenMatrix Klocal = mCachedLocalStiffnessMatrix[i];

			// compute the deformation gradient and its polar decomposition
			Matrix3R R = mCachedRotationMatrix[i];

			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (size_t j = 0; j < GetNumLocalNodes(); j++)
			{
				size_t jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, j)];
				if (jGlobal < mNumBCs)
					continue;
				int jOffset = (jGlobal - mNumBCs) * NUM_POS_COMPONENTS;
				Vector3R f;
				for (size_t k = 0; k < GetNumLocalNodes(); k++)
				{
					size_t kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, k)];
					// construct the 3x3 block
					Matrix3R Kblock;
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							Kblock(x, y) = Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
						}
					}
					f += Kblock * (!R * GetDeformedPosition(mTetMesh->GetGlobalIndex(i, k)) - mReferencePositions[kGlobal]);
				}

				// rotate the force and add it to the global one
				f = R * f;
				for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
				{
					elasticForce[jOffset + x] += f[x];
				}
			}
		}
	}

	void FemPhysicsCorotationalElasticity::ComputeRotations()
	{
		mCachedRotationMatrix.resize(GetNumElements());
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			// compute the deformation gradient and its polar decomposition
			Matrix3R F, R, U;
			ComputeDeformationGradient(i, F);
			ComputePolarDecomposition(F, R, U);

			mCachedRotationMatrix[i] = R;
		}
	}

	// assemble the global stiffness matrix for corotational elasticity (using BB shape functions)
	void FemPhysicsCorotationalElasticity::AssembleStiffnessMatrix()
	{
		PROFILE_SCOPE("CR assemble");
		size_t numNodes = GetNumFreeNodes();
		size_t numDofs = numNodes * NUM_POS_COMPONENTS;
		if (mDenseK.rows() != numDofs)
			mDenseK.resize(numDofs, numDofs);
		mDenseK.setZero();
		EigenMatrix Rlocal;
		Rlocal.resize(GetNumLocalNodes() * NUM_POS_COMPONENTS, GetNumLocalNodes() * NUM_POS_COMPONENTS);
		mCachedRotationMatrix.resize(GetNumElements());
		for (size_t i = 0; i < (size_t)mTetMesh->GetNumElements(); i++)
		{
			EigenMatrix Klocal = mCachedLocalStiffnessMatrix[i];
			Matrix3R R = mCachedRotationMatrix[i];

			// convert the 3x3 rotation matrix to Eigen
			EigenMatrix Reig = Matrix3ToEigen<EigenMatrix>(R);

			// build the element matrix with diagonal blocks equal to R
			Rlocal.setZero();
			for (size_t j = 0; j < GetNumLocalNodes(); j++)
			{
				Rlocal.block<NUM_POS_COMPONENTS, NUM_POS_COMPONENTS>(j * NUM_POS_COMPONENTS, j * NUM_POS_COMPONENTS) = Reig;
			}

			// compute the element locally rotated stiffness matrix
			EigenMatrix Klocal_rot = Rlocal * Klocal * Rlocal.transpose();

			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (size_t j = 0; j < GetNumLocalNodes(); j++)
			{
				size_t jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, j)];
				if (jGlobal < mNumBCs)
					continue;
				for (size_t k = 0; k < GetNumLocalNodes(); k++)
				{
					size_t kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, k)];
					if (kGlobal < mNumBCs)
						continue;
					int jOffset = (jGlobal - mNumBCs) * NUM_POS_COMPONENTS;
					int kOffset = (kGlobal - mNumBCs) * NUM_POS_COMPONENTS;

					// add the the whole 3x3 block to the block matrix
					for (uint32 x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (uint32 y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							mDenseK.coeffRef(jOffset + x, kOffset + y) += Klocal_rot(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
						}
					}
				}
			}
		}

		mStiffnessMatrix = mDenseK.sparseView();
	}

	void FemPhysicsCorotationalElasticity::CacheLocalStiffnessMatrices()
	{
		mCachedLocalStiffnessMatrix.resize(GetNumElements());
		for (size_t i = 0; i < (size_t)mTetMesh->GetNumElements(); i++)
		{
			EigenMatrix Klocal;
			const real mu = GetShearModulus();
			const real lambda = GetLameFirstParam();
			ComputeLocalStiffnessMatrixBB(i, mu, lambda, Klocal);
			mCachedLocalStiffnessMatrix[i] = Klocal;
		}
	}


} // namespace FEM_SYSTEM
