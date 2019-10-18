#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "FemDataStructures.h"

#pragma warning( disable : 4305)

namespace FEM_SYSTEM
{
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
		void InitSparse(SparseMatrix& matrix, LinearSolverType type);
		EigenVector Solve(const EigenVector& rhs);

	private:
		LinearSolverType mType;
		bool mSparse = false;

		// matrix decompositions
		Eigen::FullPivLU<EigenMatrix> mFullPivLU;
		Eigen::PartialPivLU<EigenMatrix> mLU;
		Eigen::LLT<EigenMatrix> mLLT;
		Eigen::LDLT<EigenMatrix> mLDLT;

		// sparse matrix decompositions
		Eigen::SimplicialLDLT<SparseMatrix> mSimplicialLDLT;
		Eigen::SparseLU<SparseMatrix> mSparseLU;
#if !defined(_DEBUG) && defined(USE_MKL)
		Eigen::PardisoLDLT<SparseMatrix> mPardisoLU;
		Eigen::PardisoLU<SparseMatrix> mPardisoLDLT;
#endif

		// iterative solvers
		Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<real>> mCG;
		Eigen::GMRES<SparseMatrix, Eigen::IncompleteLUT<real>> mGMRES;
		Eigen::MINRES<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<real>> mMINRES;
	};

} // namespace FEM_SYSTEM

#endif // LINEAR_SOLVER_H