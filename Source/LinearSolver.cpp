#include "LinearSolver.h"

namespace FEM_SYSTEM
{
	void LinearSolver::Init(const EigenMatrix& matrix, LinearSolverType type, bool sparse)
	{
		mType = type;
		mSparse = sparse;

		// compute decomposition
		if (mType == LST_LLT)
		{
			mLLT = matrix.llt();
		}
		else if (mType == LST_LU)
		{
			if (sparse)
				mSparseLU.compute(matrix.sparseView());
			else
				mLU = matrix.lu();
		}
		else if (mType == LST_LDLT)
		{
			if (mSparse)
				mSimplicialLDLT.compute(matrix.sparseView());
			else
				mLDLT = matrix.ldlt();
		}
		else if (mType == LST_FULL_LU)
		{
			mFullPivLU = matrix.fullPivLu();
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

	void LinearSolver::InitSparse(SparseMatrix& matrix, LinearSolverType type)
	{
		mType = type;
		mSparse = true;

		// compute decomposition
		if (mType == LST_LU)
		{
			matrix.makeCompressed();
			mSparseLU.compute(matrix);
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
			sol = mCG.solve(rhs);
		}
		else if (mType == LST_MINRES)
		{
			sol = mMINRES.solve(rhs);
		}
		else if (mType == LST_GMRES)
		{
			sol = mGMRES.solve(rhs);
		}

		return sol;
	}

}