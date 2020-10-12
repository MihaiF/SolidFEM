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

#include "SaddlePointProblem.h"
#include "FemDataStructures.h"

namespace FEM_SYSTEM
{
	EigenVector SaddlePointProblem::ComputeRhs(const EigenVector& solution)
	{
		uint32 numNodes = mFPM.GetNumFreeNodes();
		uint32 numPNodes = mFPM.GetNumPressureNodes();
		uint32 numContacts = mFPM.GetNumContacts();
		uint32 numConstraints = mFPM.mCondenseContact ? numPNodes : numPNodes + numContacts;
		uint32 numDofs = numNodes * NUM_POS_COMPONENTS;
		uint32 size = numDofs + numConstraints;

		real scale = mFPM.mConfig.mConstraintScale; // scaling of volume constraints

		// update deformed positions
		Vector3Array u = GetStdVector(solution.head(numDofs));
		for (uint32 i = 0; i < mFPM.GetNumFreeNodes(); i++)
		{
			mFPM.mDeformedPositions[i] = u[i];
		}

		// update pressure
		auto p = solution.segment(numDofs, numPNodes);
		for (uint32 i = 0; i < numPNodes; i++)
		{
			mFPM.mPressures[i] = p[i];
		}

		if (!mFPM.mCondenseContact)
		{
			// update contact multipliers
			mFPM.mContactMultipliers = solution.tail(numContacts);
		}

		// assemble Jacobian matrix
		mFPM.AssembleJacobianMatrix(mFPM.mVolJacobianMatrix, false, true);
#ifdef CHECK_JACOBIAN
		std::cout << "J analytic" << std::endl << mFPM.mVolJacobianMatrix << std::endl;
		EigenMatrix J = mFPM.mVolJacobianMatrix;
		mFPM.AssembleJacobianMatrixFD(J);
		std::cout << "J approx" << std::endl << J << std::endl;
#endif

		// assemble right hand side
		EigenVector rhs(size);
		rhs.setZero();

		rhs.head(numDofs) = mFPM.ComputePosRhs();

		// start with the pressure components
		EigenVector errors(numPNodes); // vector used to accumulate volume errors
		real totalErr = mFPM.GetCurrentVolumeErrors(errors);
		rhs.segment(numDofs, numPNodes) = -scale * errors + scale * mFPM.mVolComplianceMatrix * p;

		// continue with contacts
		if (!mFPM.mCondenseContact && numContacts != 0)
			rhs.tail(numContacts) = -mFPM.mContactDepth;

		return rhs;
	}

	real SaddlePointProblem::MeritResidual(const EigenVector& rhs) const
	{
		return 0.5 * rhs.squaredNorm();
	}

	void SaddlePointProblem::ComputeSystemMatrix(const EigenVector& s, const EigenVector& y, EigenMatrix& K)
	{
		if (mSolverType != SPP_ST_DIRECT)
			return;

		uint32 numNodes = mFPM.GetNumFreeNodes();
		uint32 numPNodes = mFPM.GetNumPressureNodes();
		uint32 numContacts = mFPM.GetNumContacts();
		uint32 numConstraints = mFPM.mCondenseContact ? numPNodes : numPNodes + numContacts;
		uint32 numDofs = numNodes * NUM_POS_COMPONENTS;
		uint32 size = numDofs + numConstraints;

		real scale = mFPM.mConfig.mConstraintScale; // scaling of volume constraints

		// the system matrix
		K.resize(size, size);

		auto Kxx = K.block(0, 0, numDofs, numDofs);
		mFPM.ComputeStiffnessMatrix(Kxx);
		K.block(0, numDofs, numDofs, numConstraints) = scale * mFPM.mVolJacobianMatrix.transpose();
		K.block(numDofs, 0, numConstraints, numDofs) = scale * mFPM.mVolJacobianMatrix;
		K.block(numDofs, numDofs, numPNodes, numPNodes) = -scale * mFPM.mVolComplianceMatrix;
	}

	void SaddlePointProblem::ComputeSystemMatrix(const EigenVector& s, const EigenVector& y, SparseMatrix& K)
	{
		if (mSolverType != SPP_ST_DIRECT)
			return;

		uint32 numNodes = mFPM.GetNumFreeNodes();
		uint32 numPNodes = mFPM.GetNumPressureNodes();
		uint32 numContacts = mFPM.GetNumContacts();
		uint32 numConstraints = mFPM.mCondenseContact ? numPNodes : numPNodes + numContacts;
		uint32 numDofs = numNodes * NUM_POS_COMPONENTS;
		uint32 size = numDofs + numConstraints;

		real scale = mFPM.mConfig.mConstraintScale; // scaling of volume constraints

		// the system matrix
		K = SparseMatrix(size, size);

		std::vector<Eigen::Triplet<real>> triplets;
		SparseMatrix Kxx;
		mFPM.ComputeStiffnessMatrix(Kxx);
		for (int k = 0; k < Kxx.outerSize(); ++k)
		{
			for (SparseMatrix::InnerIterator it(Kxx, k); it; ++it)
			{
				triplets.push_back(Eigen::Triplet<real>((int)it.row(), (int)it.col(), it.value()));
			}
		}
		for (int k = 0; k < mFPM.mVolJacobianMatrix.outerSize(); ++k)
		{
			for (SparseMatrix::InnerIterator it(mFPM.mVolJacobianMatrix, k); it; ++it)
			{
				triplets.push_back(Eigen::Triplet<real>((int)(it.row() + numDofs), (int)it.col(), scale * it.value()));
				triplets.push_back(Eigen::Triplet<real>((int)it.col(), (int)(it.row() + numDofs), scale * it.value()));
			}
		}
		for (int k = 0; k < mFPM.mVolComplianceMatrix.outerSize(); ++k)
		{
			for (SparseMatrix::InnerIterator it(mFPM.mVolComplianceMatrix, k); it; ++it)
			{
				triplets.push_back(Eigen::Triplet<real>((int)(it.row() + numDofs), (int)(it.col() + numDofs), -scale * it.value()));
			}
		}
		K.setFromTriplets(triplets.begin(), triplets.end());
	}

	EigenVector SaddlePointProblem::SolveLinearSystem(SparseMatrix& A, const EigenVector& rhs, const EigenVector& s, const EigenVector& y)
	{
		if (mSolverType == SPP_ST_DIRECT)
		{
			LinearSolver solver;
			solver.Init(A, LST_LDLT_PARDISO);
			return solver.Solve(rhs);
		}
		return EigenVector();
	}

	EigenVector SaddlePointProblem::SolveLinearSystem(EigenMatrix& A, const EigenVector& rhs, const EigenVector& s, const EigenVector& y)
	{
		// TODO: no real need for K here

		uint32 numNodes = mFPM.GetNumFreeNodes();
		uint32 numPNodes = mFPM.GetNumPressureNodes();
		uint32 numContacts = mFPM.GetNumContacts();
		uint32 numConstraints = mFPM.mCondenseContact ? numPNodes : numPNodes + numContacts;
		uint32 numDofs = numNodes * NUM_POS_COMPONENTS;
		uint32 size = numDofs + numConstraints;
		real invHSqr = 1.f / (mFPM.mTimeStep * mFPM.mTimeStep);

		if (mSolverType == SPP_ST_DIRECT)
		{
			LinearSolver solver;
			solver.Init(A, LST_LDLT_PARDISO);
			return solver.Solve(rhs);
		}
		else if (mSolverType == SPP_ST_SCHUR_COMPLIANCE)
		{
			// static condensation
			EigenMatrix K;
			mFPM.ComputeStiffnessMatrix(K);
			K += mFPM.mVolJacobianMatrix.transpose() * mCinv * mFPM.mVolJacobianMatrix;
			if (mFPM.mSimType == ST_IMPLICIT)
				K += invHSqr * mFPM.mMassMatrix;

			LinearSolver solver;
			solver.Init(K, LST_LU_PARDISO);

			EigenVector solution(size);
			solution.head(numDofs) = solver.Solve(rhs.head(numDofs) + mFPM.mVolJacobianMatrix.transpose() * mCinv * rhs.segment(numDofs, numPNodes));
			solution.segment(numDofs, numPNodes) = mCinv * (mFPM.mVolJacobianMatrix * solution.head(numDofs) - rhs.segment(numDofs, numPNodes)); // back-solve
			return solution;
		}
		else if (mSolverType == SPP_ST_SCHUR_MASS)
		{
			{
				EigenMatrix K = mFPM.mDeviatoricStiffnessMatrix + mFPM.mGeometricStiffnessMatrix;
				if (mFPM.mSimType == ST_IMPLICIT)
					K += invHSqr * mFPM.mMassMatrix;
				mKinv = K.inverse(); // compute inverse
			}

			auto JKinv = mFPM.mVolJacobianMatrix * mKinv;
			EigenMatrix S = JKinv * mFPM.mVolJacobianMatrix.transpose() + mFPM.mVolComplianceMatrix;

			EigenVector solution(size);

			// solve reduced system
			LinearSolver solver;
			solver.Init(S, LST_LU_PARDISO);
			solution.tail(numPNodes) = solver.Solve(-rhs.tail(numPNodes) + JKinv * rhs.head(numDofs));

			// back-solve
			solution.head(numDofs) = mKinv * (rhs.head(numDofs) - mFPM.mVolJacobianMatrix.transpose() * solution.tail(numPNodes));

			return solution;
		}
		else if (mSolverType == SPP_ST_ITERATIVE)
		{
			std::vector<real> stdFullSol(size);
			std::vector<real> stdRhs(rhs.data(), rhs.data() + size);
			int maxIterations = 100;
			real tolerance = 0.1; // this is actually the forcing term from inexact Newton (?)

			LinearSolver solver;
			//solver.SetTolerance(tolerance);
			solver.Init(A, LST_CG);
			auto solution = solver.Solve(rhs);
			Printf("Krylov iterations: %d\n", solver.GetIterations());
			return solution;
		}

		ASSERT(false);
		return EigenVector();
	}
}