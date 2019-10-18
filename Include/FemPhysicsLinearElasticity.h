#ifndef FEM_PHYSICS_LINEAR_EXPLICIT_H
#define FEM_PHYSICS_LINEAR_EXPLICIT_H

#include "FemPhysicsLinear.h"
#include "LinearSolver.h"

namespace FEM_SYSTEM
{
	class FemPhysicsLinearElasticity : public FemPhysicsLinear
	{
	public:
		FemPhysicsLinearElasticity(std::vector<Tetrahedron>& tetrahedra,
			std::vector<Node>& nodes, const FemConfig& config);
		void Step(real dt) override;
		void SolveEquilibrium(float) override;

	private:
		void StepExplicit(real h);
		void StepImplicit(real h);
		void AssembleStiffnessMatrix();

	protected:
		// linear elasticity helpers
		void ComputeLocalStiffnessMatrix(uint32 i, EigenMatrix& Klocal);
		void ComputeLocalStiffnessMatrixBB1(uint32 i, EigenMatrix& Klocal);
		void ComputeDyadicMatrixBB2(uint32 i, uint32 j, Vector3R y[4], Matrix3R& L);
		void ComputeLocalStiffnessMatrixBB2(uint32 i, EigenMatrix& Klocal);

	protected:
		EigenMatrix mStiffnessMatrix; // stiffness matrix
		Eigen::SparseLU<SparseMatrix> mSparseLU; // the LU decomposition of M
		
		// for implicit integration
		bool mHaveSysMatrix = false;
		LinearSolver mSolver;

		friend class FemTester;
	};

} // namespace FEM_SYSTEM

#endif // !FEM_PHYSICS_LINEAR_EXPLICIT_H
