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

#include "LinearSolver.h"
#include <Engine/Utils.h>

//namespace FEM_SYSTEM
//{
	void LinearSolver::Init(const EigenMatrix& matrix, LinearSolverType type, bool sparse)
	{
		mType = type;
		mSparse = sparse;

#if defined(_DEBUG) || !defined(USE_MKL)
		// fall back on Eigen if MKL/Pardiso support absent
		if (mType == LST_LDLT_PARDISO)
			mType = LST_LDLT;
		if (mType == LST_LU_PARDISO)
			mType = LST_LU;
#endif

		// compute decomposition
		if (mType == LST_LLT)
		{
			mLLT = matrix.llt();
			if (mLLT.info() != Eigen::Success)
				Printf("LLT failed\n");
		}
		else if (mType == LST_LU)
		{
			if (sparse)
				mSparseLU.compute(matrix.sparseView());
			else
			{
				mLU = matrix.lu();
				Printf("det: %g\n", mLU.determinant());
			}
		}
		else if (mType == LST_LDLT)
		{
			if (mSparse)
				mSimplicialLDLT.compute(matrix.sparseView());
			else
			{
				mLDLT = matrix.ldlt();
				Printf("PD: %d\n", mLDLT.isPositive());
			}
		}
		else if (mType == LST_FULL_LU)
		{
			mFullPivLU = matrix.fullPivLu();
			Printf("det: %g\n", mFullPivLU.determinant());
			Printf("rank: %d / %d\n", mFullPivLU.rank(), matrix.rows());
			if (!mFullPivLU.isInvertible())
				Printf("Matrix is singular\n");
		}
#if !defined(_DEBUG) && defined(USE_MKL)
		else if (mType == LST_LDLT_PARDISO)
		{
			mPardisoLDLT.compute(matrix.sparseView());
		}
		else if (mType == LST_LU_PARDISO)
		{
			mPardisoLU.compute(matrix.sparseView());
		}
#endif
		else if (mType == LST_CG)
		{
			mCG.compute(matrix.sparseView());
		}
		else if (mType == LST_MINRES)
		{
			mMINRES.compute(matrix.sparseView());
		}
		else if (mType == LST_GMRES)
		{
			mGMRES.compute(matrix.sparseView());
		}
	}

	void LinearSolver::Init(SparseMatrix& matrix, LinearSolverType type)
	{
		mType = type;
		mSparse = true;

#if defined(_DEBUG) || !defined(USE_MKL)
		// fall back on Eigen if MKL/Pardiso support absent
		if (mType == LST_LDLT_PARDISO)
			mType = LST_LDLT;
		if (mType == LST_LU_PARDISO)
			mType = LST_LU;
#endif

		// compute decomposition
		if (mType == LST_LU || mType == LST_FULL_LU)
		{
			matrix.makeCompressed();
			mSparseLU.compute(matrix);
			Printf("det: %g\n", mSparseLU.determinant());
		}
		else if (mType == LST_LLT)
		{
			mSimplicialLLT.compute(matrix);
			if (mSimplicialLLT.info() != Eigen::Success)
				Printf("LLT failed\n");
		}
		else if (mType == LST_LDLT)
		{
			mSimplicialLDLT.compute(matrix);
		}
#if !defined(_DEBUG) && defined(USE_MKL)
		else if (mType == LST_LDLT_PARDISO)
		{
			mPardisoLDLT.compute(matrix);
		}
		else if (mType == LST_LU_PARDISO)
		{
			mPardisoLU.compute(matrix);
		}
#endif
		else if (mType == LST_CG)
		{
			mCG.compute(matrix);
		}
		else if (mType == LST_MINRES)
		{
			mMINRES.compute(matrix);
		}
		else if (mType == LST_GMRES)
		{
			mGMRES.compute(matrix);
		}
	}

	EigenVector LinearSolver::Solve(const EigenVector& rhs)
	{
		EigenVector sol;
		if (mType == LST_LLT)
		{
			sol = mLLT.solve(rhs);
		}
		else if (mType == LST_LU)
		{
			if (mSparse)
				sol = mSparseLU.solve(rhs);
			else
				sol = mLU.solve(rhs);
		}
		else if (mType == LST_LDLT)
		{
			if (mSparse)
				sol = mSimplicialLDLT.solve(rhs);
			else
				sol = mLDLT.solve(rhs);
		}
		else if (mType == LST_FULL_LU)
		{
			sol = mFullPivLU.solve(rhs);
		}
#if !defined(_DEBUG) && defined(USE_MKL)
		else if (mType == LST_LDLT_PARDISO)
		{
			sol = mPardisoLDLT.solve(rhs);
		}
		else if (mType == LST_LU_PARDISO)
		{
			sol = mPardisoLU.solve(rhs);
		}
#endif
		else if (mType == LST_CG)
		{
			mCG.setTolerance(mTolerance);
			if (mMaxIterations > 0)
				mCG.setMaxIterations(mMaxIterations);
			sol = mCG.solve(rhs);
			mIterations = (int)mCG.iterations();
			//Printf("CG statistics: iterations=%d / %d; error=%g / %g\n",
			//	mCG.iterations(), mCG.maxIterations(), mCG.error(), mCG.tolerance());
		}
		else if (mType == LST_MINRES)
		{
			mMINRES.setTolerance(mTolerance);
			sol = mMINRES.solve(rhs);
			mIterations = (int)mMINRES.iterations();
		}
		else if (mType == LST_GMRES)
		{
			//mGMRES.setTolerance(mTolerance);
			sol = mGMRES.solve(rhs);
			mIterations = (int)mGMRES.iterations();
		}

		return sol;
	}

//} //namespace