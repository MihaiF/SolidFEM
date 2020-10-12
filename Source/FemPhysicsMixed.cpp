/*
BSD 3-Clause License

Copyright (c) 2020, Mihai Francu
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

#include "FemPhysicsMixed.h"
#include <Engine/Profiler.h>
#include "ElasticEnergy.h"
#include "NewtonSolver.h"
#include "MixedProblems.h"
#include "SaddlePointProblem.h"

namespace FEM_SYSTEM
{

	FemPhysicsMixed::FemPhysicsMixed(std::vector<Tet>& tetrahedra, std::vector<Node>& allNodes, const FemConfig& config) 
		: FemPhysicsLinearIncompressible(tetrahedra, allNodes, config)
	{
		std::fill(mPressures.begin(), mPressures.end(), 0); // zero pressures
	}

	void FemPhysicsMixed::Step(real dt)
	{
		mPreviousPositions = mDeformedPositions;

		if (mSimType == ST_STATIC)
		{
			if (mForceFraction == 0)
			{
				mForceFraction = 1;
				SolveSaddlePoint();
			}
		}
		else if (mSimType == ST_QUASI_STATIC)
		{
			if (mForceFraction < 1)
			{
				mForceFraction = std::min(real(1), mForceFraction + mForceStep);
				if (mConfig.mSolver == NST_MIXED_NEWTON || mConfig.mSolver == NST_NEWTON)
					SolveSaddlePoint();
				else if (mConfig.mSolver == NST_NEWTON_LS)
					SolveSaddlePointLS();

				EigenVector errors(GetNumPressureNodes());
				GetCurrentVolumeErrors(errors, true);
			}
		}
		else if (mSimType == ST_IMPLICIT)
		{
			mForceFraction = 1;
			real h = dt / mNumSteps;
			for (int i = 0; i < mNumSteps; i++)
			{
				mTimeStep = h;
				SolveSaddlePoint();
				HandleCollisions(h);
			}
		}
	}

	void FemPhysicsMixed::SolveSaddlePoint()
	{
		MEASURE_TIME("Mixed FEM solve");
		PROFILE_SCOPE("Mixed");

		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumPressureNodes();
		uint32 numContacts = GetNumContacts();
		uint32 numConstraints = mCondenseContact ? numPNodes : numPNodes + numContacts;
		uint32 numDofs = numNodes * NUM_POS_COMPONENTS;
		uint32 size = numDofs + numConstraints;

		AssembleDynamicContributions();

		NewtonSolver<SaddlePointProblem, SparseMatrix> solver;
		solver.mNumIterations = mOuterIterations;
		solver.mVerbose = mVerbose ? VL_DETAILED : VL_NONE;
		solver.mResidualThreshold = mAbsNewtonResidualThreshold;
		solver.mSolverType = LST_LDLT_PARDISO;
		SaddlePointProblem problem(*this, SaddlePointProblem::SPP_ST_DIRECT);

		// prepare the initial guess with the current configuration
		EigenVector solution(size);
		solution.setZero(); // initialize Lagrange multipliers with zero
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			const Vector3R& p = mDeformedPositions[i];
			solution(i * NUM_POS_COMPONENTS) = p.x;
			solution(i * NUM_POS_COMPONENTS + 1) = p.y;
			solution(i * NUM_POS_COMPONENTS + 2) = p.z;
		}
		// warm-starting
		for (uint32 i = 0; i < numPNodes; i++)
		{
			solution(numDofs + i) = mPressures[i];
		}

		solver.Solve(problem, size, solution);

		// update deformed positions
		Vector3Array u = GetStdVector(solution.head(numDofs));
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			mDeformedPositions[i] = u[i];
		}

		// update pressure
		auto p = solution.tail(numPNodes);
		for (uint32 i = 0; i < numPNodes; i++)
		{
			mPressures[i] = p[i];
		}

		CheckForInversion(true);

		// compute new velocities
		if (mSimType == ST_IMPLICIT)
		{
			real invH = 1.f / mTimeStep;
			for (uint32 i = 0; i < GetNumFreeNodes(); i++)
			{
				mVelocities[i] = invH * (mDeformedPositions[i] - mPreviousPositions[i]);
			}
		}
	}

	void FemPhysicsMixed::SolveSaddlePointLS()
	{
		MEASURE_TIME("Mixed FEM solve");
		PROFILE_SCOPE("Mixed");

		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumPressureNodes();
		uint32 numContacts = GetNumContacts();
		uint32 numConstraints = mCondenseContact ? numPNodes : numPNodes + numContacts;
		uint32 numDofs = numNodes * NUM_POS_COMPONENTS;
		uint32 size = numDofs + numConstraints;

		AssembleDynamicContributions();

		NewtonSolverBackTrack<SaddlePointProblem, SparseMatrix> solver;
		solver.mNumIterations = mOuterIterations;
		solver.mVerbose = mVerbose ? VL_MINIMUM : VL_NONE;
		solver.mResidualThreshold = mAbsNewtonResidualThreshold;
		solver.mUseProblemSolver = true;
		SaddlePointProblem problem(*this, SaddlePointProblem::SPP_ST_DIRECT);

		// prepare the initial guess with the current configuration
		EigenVector solution(size);
		solution.setZero(); // initialize Lagrange multipliers with zero
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			const Vector3R& p = mDeformedPositions[i];
			solution(i * NUM_POS_COMPONENTS) = p.x;
			solution(i * NUM_POS_COMPONENTS + 1) = p.y;
			solution(i * NUM_POS_COMPONENTS + 2) = p.z;
		}
		// warm-starting
		for (uint32 i = 0; i < numPNodes; i++)
		{
			solution(numDofs + i) = mPressures[i];
		}

		solver.Solve(problem, size, solution);

		// update deformed positions
		Vector3Array u = GetStdVector(solution.head(numDofs));
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			mDeformedPositions[i] = u[i];
		}

		// update pressure
		auto p = solution.tail(numPNodes);
		for (uint32 i = 0; i < numPNodes; i++)
		{
			mPressures[i] = p[i];
		}

		CheckForInversion(true);

		// compute new velocities
		if (mSimType == ST_IMPLICIT)
		{
			real invH = 1.f / mTimeStep;
			for (uint32 i = 0; i < GetNumFreeNodes(); i++)
			{
				mVelocities[i] = invH * (mDeformedPositions[i] - mPreviousPositions[i]);
			}
		}
	}

	real FemPhysicsMixed::GetCurrentVolumeErrors(EigenVector& errors, bool verbose)
	{
		real totalVol = 0;
		mTotalInitialVol = 0;
		real totalErr = 0;
		real factor = mPressureOrder == 1 ? 0.25f : 1;
		errors.setZero();
		for (uint32 e = 0; e < GetNumElements(); e++)
		{
			const Vector3R& x0 = GetDeformedPosition(mTetMesh->GetGlobalIndex(e, 0));
			const Vector3R& x1 = GetDeformedPosition(mTetMesh->GetGlobalIndex(e, 1));
			const Vector3R& x2 = GetDeformedPosition(mTetMesh->GetGlobalIndex(e, 2));
			const Vector3R& x3 = GetDeformedPosition(mTetMesh->GetGlobalIndex(e, 3));
			Vector3R d1 = x1 - x0;
			Vector3R d2 = x2 - x0;
			Vector3R d3 = x3 - x0;
			Matrix3R mat(d1, d2, d3); // this is the spatial shape matrix Ds [Sifakis][Teran]
			real vol = (mat.Determinant()) / 6.f; // volume of the tet
			mTotalInitialVol += mElementVolumes[e];
			real err;
			if (mLogConstraint)
			{
				real J = vol / mElementVolumes[e];
				err = mElementVolumes[e] * log(J);
			}
			else
				err = vol - mElementVolumes[e];
			totalVol += vol;
			totalErr += err;

			for (uint32 i = 0; i < GetNumLocalPressureNodes(); i++)
			{
				uint32 globalI = GetPressureGlobalIndex(e, i);
				errors[globalI] += factor * err;
			}
		}
		if (verbose)
		{
			Printf("vol err: %.4f%%\n", abs(totalErr) / mTotalInitialVol * 100);
		}
		return totalErr / mTotalInitialVol;
	}

	void FemPhysicsMixed::AssembleGeometricStiffnessMatrixFD(const EigenVector& p, EigenMatrix& Kgeom)
	{
		uint32 numNodes = GetNumFreeNodes();
		uint32 numDofs = numNodes * 3;
		Kgeom.resize(numDofs, numDofs);
		Kgeom.setZero();

		// current pressure forces
		auto fp0 = mVolJacobianMatrix.transpose() * p;

		auto J = mVolJacobianMatrix;
		real eps = 1e-6;
		for (uint32 j = 0; j < numDofs; j++)
		{
			mDeformedPositions[j / 3][j % 3] += eps;
			AssembleJacobianMatrix(J, false, true);
			auto fp = J.transpose() * p;
			auto dfp = (1.0 / eps) * (fp - fp0);
			for (uint32 i = 0; i < numDofs; i++)
			{
				Kgeom(i, j) = dfp(i);
			}
			mDeformedPositions[j / 3][j % 3] -= eps;
		}
	}

	void FemPhysicsMixed::AssembleGeometricStiffnessMatrix(const EigenVector& p, SparseMatrix& Kgeom) const
	{
		PROFILE_SCOPE("Geom stiff");

		size_t numNodes = GetNumFreeNodes();
		size_t numDofs = numNodes * 3;
		Kgeom.resize(numDofs, numDofs);
		Kgeom.setZero();
		std::vector<Eigen::Triplet<real>> triplets;
		
		uint32 numPressureLocalNodes = GetNumLocalPressureNodes();
		real factor = mPressureOrder == 1 ? 0.25f : 1.f;

		// go through all linear elements (tetrahedra)
		Matrix3R Z = Matrix3R::Zero();
		Eigen::Matrix<real, 12, 12> Klocal;
		for (int e = 0; e < (int)GetNumElements(); e++)
		{
			// compute local geometric stiffness matrix
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

			// prepare the building blocks
			Matrix3R K1 = Matrix3R::Skew(d1);
			Matrix3R K2 = Matrix3R::Skew(d2);
			Matrix3R K3 = Matrix3R::Skew(d3);

			// assemble the block matrix
			Matrix3R K[4][4];
			K[0][0] = Z; K[0][1] = K2 - K3; K[0][2] = K3 - K1; K[0][3] = K1 - K2;
			K[1][0] = K3 - K2; K[1][1] = Z; K[1][2] = -K3; K[1][3] = K2;
			K[2][0] = -K3 + K1; K[2][1] = K3; K[2][2] = Z; K[2][3] = -K1;
			K[3][0] = K2 - K1; K[3][1] = -K2; K[3][2] = K1; K[3][3] = Z;

			// convert it to an EigenMatrix
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					for (int x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (int y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							Klocal(i * NUM_POS_COMPONENTS + x, j * NUM_POS_COMPONENTS + y) = K[i][j](x, y);
						}
					}
				}
			}

			real factor1 = factor / 6;
			if (mLogConstraint)
			{
				Matrix3R mat(d1, d2, d3);
				real vol = mat.Determinant() / 6;
				real J = vol / mElementVolumes[e];
				factor1 = factor1 / J;

				Vector3R v[4];
				ComputeLocalJacobian(e, v, true);

				EigenMatrix Jlocal(1, 12);
				Jlocal.setZero();
				{
					// build Jacobian matrix
					for (uint32 k = 0; k < GetNumLocalNodes(); k++)
					{
						for (int l = 0; l < NUM_POS_COMPONENTS; l++)
						{
							Jlocal(0, k * NUM_POS_COMPONENTS + l) = v[k][l];
						}
					}
				}
				Klocal += -(1 / vol) * Jlocal.transpose() * Jlocal;
			}

			for (uint32 i = 0; i < numPressureLocalNodes; i++)
			{
				uint32 globalI = GetPressureGlobalIndex(e, i);
				// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
				for (uint32 j = 0; j < GetNumLocalNodes(); j++)
				{
					uint32 jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(e, (int)j)];
					for (uint32 k = 0; k < GetNumLocalNodes(); k++)
					{
						uint32 kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(e, (int)k)];
						if (jGlobal < mNumBCs || kGlobal < mNumBCs)
							continue;
						uint32 jOffset = (jGlobal - mNumBCs) * NUM_POS_COMPONENTS;
						uint32 kOffset = (kGlobal - mNumBCs) * NUM_POS_COMPONENTS;

						// add the the whole 3x3 block to the block matrix
						for (uint32 x = 0; x < NUM_POS_COMPONENTS; x++)
						{
							for (uint32 y = 0; y < NUM_POS_COMPONENTS; y++)
							{
								real val = factor1 * p[globalI] * Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
								triplets.push_back(Eigen::Triplet<real>(jOffset + x, kOffset + y, val));
							}
						}
					}
				}
			}
		}
		Kgeom.setFromTriplets(triplets.begin(), triplets.end());
	}

	void FemPhysicsMixed::ComputeGradients(std::vector<Vector3R>& r, const EigenVector& dx, int material)
	{
		// compute negative gradient r = f(x) + M * ( v / h - dx / h^2 )
		for (size_t i = 0; i < r.size(); i++)
		{
			r[i].SetZero();
		}
		ElasticEnergy::ComputeForces(this, r, material);
		if (mSimType == ST_IMPLICIT)
		{
			const real invH = 1.f / mTimeStep;
			const real invHSqr = 1.f / (mTimeStep * mTimeStep);
			EigenVector inertia = mMassMatrix * (invH * GetEigenVector(mVelocities) - invHSqr * dx);
			auto stdInertia = GetStdVector(inertia);
			for (size_t i = 0; i < GetNumFreeNodes(); i++)
			{
				r[i + mNumBCs] += stdInertia[i];
			}
		}
	}

	template<class MATRIX>
	void FemPhysicsMixed::ComputeStiffnessMatrix(MATRIX& K, bool update)
	{
		real scale = mConfig.mConstraintScale;
		real invHSqr = 1.f / (mTimeStep * mTimeStep);
		
		if (update)
		{
			Eigen::Map<EigenVector> p((real*)&mPressures[0], mPressures.size(), 1);
			AssembleGeometricStiffnessMatrix(p, mGeometricStiffnessMatrix);
			// contacts and Dirichlet BCs have no geometric stiffness because they are linear constraints, i.e. the Hessian vanishes
			ElasticEnergy::AssembleStiffnessMatrix(this, mDeviatoricStiffnessMatrix);
		}

		K = mDeviatoricStiffnessMatrix + scale * mGeometricStiffnessMatrix;
		if (mSimType == ST_IMPLICIT)
			K += invHSqr * mMassMatrix;
		if (GetNumContacts() != 0 && mCondenseContact)
		{
			K += mContactStiffness * mContactJacobian.transpose() * mContactJacobian;
		}
		if (!mDirichletIndices.empty())
		{
			K += mDirichletStiffness * mDirichletJacobian.transpose() * mDirichletJacobian;
		}
	}

	// explicit instantiations
	template void FemPhysicsMixed::ComputeStiffnessMatrix<EigenMatrix>(EigenMatrix& K, bool update);
	template void FemPhysicsMixed::ComputeStiffnessMatrix<SparseMatrix>(SparseMatrix& K, bool update);
	template void FemPhysicsMixed::ComputeStiffnessMatrix<Eigen::Block<EigenMatrix>>(Eigen::Block<EigenMatrix>& K, bool update);

	EigenVector FemPhysicsMixed::ComputeStdRhs(int material)
	{
		uint32 numDofs = GetNumFreeNodes() * NUM_POS_COMPONENTS;

		EigenVector rhs(numDofs);

		EigenVector f = mForceFraction * GetEigenVector(mBodyForces);
		rhs = f;

		ComputeTractionForces();
		if (mTractionForces.size() > 0)
			rhs += mForceFraction * GetEigenVector(mTractionForces);

		// elastic forces vector
		std::vector<Vector3R> fe(GetNumNodes());
		EigenVector disp = GetEigenVector(mDeformedPositions) - GetEigenVector(mPreviousPositions);
		ComputeGradients(fe, disp, material);
		rhs += GetEigenVector(fe, mNumBCs);

		return rhs;
	}

	EigenVector FemPhysicsMixed::ComputePosRhs()
	{
		uint32 numDofs = GetNumFreeNodes() * NUM_POS_COMPONENTS;
		uint32 numContacts = GetNumContacts();
		real scale = mConfig.mConstraintScale; // scaling of volume constraints

		EigenVector rhs = ComputeStdRhs();

		// add current pressure forces
		Eigen::Map<EigenVector> p((real*)&mPressures[0], mPressures.size(), 1);
		rhs.head(numDofs) -= scale * mVolJacobianMatrix.transpose() * p;

		// add penalty Dirichlet BC forces
		if (!mDirichletIndices.empty())
		{
			// count the constraints
			int count = 0;
			for (size_t i = 0; i < mDirichletIndices.size(); i++)
			{
				if (mDirichletAxes[i] & AXIS_X) count++;
				if (mDirichletAxes[i] & AXIS_Y) count++;
				if (mDirichletAxes[i] & AXIS_Z) count++;
			}

			EigenVector dirichletErrors(count);
			count = 0;
			for (uint32 i = 0; i < mDirichletIndices.size(); i++)
			{
				// the i'th BC				
				uint32 flags = mDirichletAxes[i];
				uint32 idx0 = mDirichletIndices[i];
				uint32 idx = idx0 - mNumBCs; // affecting node idx
				// compute the BC error
				if (flags & AXIS_X)
				{
					dirichletErrors(count++) = mDeformedPositions[idx].x - mReferencePositions[idx0].x;
				}
				if (flags & AXIS_Y)
				{
					dirichletErrors(count++) = mDeformedPositions[idx].y - mReferencePositions[idx0].y;
				}
				if (flags & AXIS_Z)
				{
					dirichletErrors(count++) = mDeformedPositions[idx].z - mReferencePositions[idx0].z;
				}
			}
			rhs.head(numDofs) -= mDirichletStiffness * mDirichletJacobian.transpose() * dirichletErrors;
		}
		return rhs;
	}

	void FemPhysicsMixed::AssembleJacobianMatrixFD(EigenMatrix& J)
	{
		uint32 numNodes = GetNumFreeNodes();
		uint32 nDof = NUM_POS_COMPONENTS * numNodes; // position DOFs
		uint32 numLocalNodes = GetNumLocalNodes();
		uint32 nLocalDofs = numLocalNodes * NUM_POS_COMPONENTS;

		uint32 numPressureLocalNodes = GetNumLocalPressureNodes();
		uint32 numPNodes = GetNumFreePressureNodes();

		J.resize(numPNodes, nDof);
		J.setZero();

		// base constraint value
		EigenVector f0(numPNodes);
		GetCurrentVolumeErrors(f0);

		real eps = 1e-6;
		EigenVector f(numPNodes);
		for (uint32 j = 0; j < nDof; j++)
		{
			mDeformedPositions[j / 3][j % 3] += eps;
			GetCurrentVolumeErrors(f);
			auto df = (1.0 / eps) * (f - f0);
			for (uint32 i = 0; i < numPNodes; i++)
			{
				J(i, j) = df(i);
			}
			mDeformedPositions[j / 3][j % 3] -= eps;
		}
	}

} // namespace FEM_SYSTEM