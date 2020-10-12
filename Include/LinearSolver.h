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

#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "FemDataStructures.h"

#pragma warning( disable : 4305)

//namespace FEM_SYSTEM
//{
	enum LinearSolverType
	{
		LST_FULL_LU,
		LST_LU,
		LST_LLT,
		LST_LDLT,
		LST_LU_PARDISO,
		LST_LDLT_PARDISO,
		LST_CG,
		LST_MINRES,
		LST_GMRES
	};

	class LinearSolver
	{
	public:
		void Init(const EigenMatrix& matrix, LinearSolverType type, bool sparse = false);
		void Init(SparseMatrix& matrix, LinearSolverType type);
		EigenVector Solve(const EigenVector& rhs);

		void SetTolerance(real tol) { mTolerance = tol; }
		void SetMaxIterations(int maxIter) { mMaxIterations = maxIter; }
		int GetIterations() const { return mIterations; }

	private:
		LinearSolverType mType;
		bool mSparse = false;
		real mTolerance = 1e-4; // for iterative solvers
		int mIterations = -1;
		int mMaxIterations = -1;

		// matrix decompositions
		Eigen::FullPivLU<EigenMatrix> mFullPivLU;
		Eigen::PartialPivLU<EigenMatrix> mLU;
		Eigen::LLT<EigenMatrix> mLLT;
		Eigen::LDLT<EigenMatrix> mLDLT;

		// sparse matrix decompositions
		Eigen::SimplicialLLT<SparseMatrix> mSimplicialLLT;
		Eigen::SimplicialLDLT<SparseMatrix> mSimplicialLDLT;
		Eigen::SparseLU<SparseMatrix> mSparseLU;
#if !defined(_DEBUG) && defined(USE_MKL)
		Eigen::PardisoLDLT<SparseMatrix> mPardisoLDLT;
		Eigen::PardisoLU<SparseMatrix> mPardisoLU;
#endif

		// iterative solvers
		Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> mCG;
		Eigen::GMRES<SparseMatrix, Eigen::IncompleteLUT<real>> mGMRES;
		Eigen::MINRES<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<real>> mMINRES;
	};

//} // namespace FEM_SYSTEM

#endif // LINEAR_SOLVER_H