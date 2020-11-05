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

#ifndef FEM_PHYSICS_LINEAR_INCOMP_H
#define FEM_PHYSICS_LINEAR_INCOMP_H

#include "FemPhysicsLinear.h"
#include "LinearSolver.h"
#include "LinearTetrahedralMesh.h"
#include <set>

#pragma warning( disable : 4305)

namespace FEM_SYSTEM
{
	class FemPhysicsLinearIncompressible : public FemPhysicsLinear
	{
	public:
		struct Config
		{
			bool mUseStaticCondensation = true;
			bool mUseKKTMatrix = true;
			bool mUseCorotational = false;
			int mPressureOrder = 1;
			bool mFullIncompressilble = false;
			real mConstraintScale = 1; // for nonlinear
			bool mLogConstraint = false;
			NonlinearSolverType mSolver = NST_MIXED_NEWTON;
		};

	public:
		FemPhysicsLinearIncompressible(std::vector<Tet>& tetrahedra,
			std::vector<Node>& nodes,
			const FemConfig& config);
		void Step(real dt) override;
		void SolveEquilibrium(float);
		void SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, const std::vector<uint32>& elemList, real pressure);

		uint32 GetNumLocalPressureNodes() const { return mPressureOrder == 0 ? 1 : 4; }
		uint32 GetNumPressureNodes() const { return mNumPressureNodes; }
		uint32 GetNumFreePressureNodes() const { return mNumPressureNodes - mNumFixedP; }
		uint32 GetPressureGlobalIndex(uint32 eidx, uint32 lidx) const { return mPressureOrder == 0 ? eidx : mPressureMesh->GetGlobalIndex(eidx, lidx); }

		void SetLameParams(real mu, real lambda) override
		{
			FemPhysicsBase::SetLameParams(mu, lambda);
			AssembleComplianceMatrix();
		}

	private:
		// explicit integration
		void StepUnconstrained(real h);
		void SolveConstraintsSchur(real h, bool implicit = false);
		void SolveConstraintsKKT(real h, bool implicit = false);

		void StepImplicitKKT(real h);
		void StepImplicitSchur(real h);
		void StepImplicitCorotationalKKT(real h);
		void StepImplicitCorotationalSchur(real h);
		
		void AssembleComplianceMatrix();

		void ComputeLocalDeviatoricStiffnessMatrix(uint32 i, EigenMatrix& Klocal);
		void ComputeLocalDeviatoricStiffnessMatrixBB2(uint32 i, EigenMatrix& Klocal);
		void ComputeLocalVolumetricStiffnessMatrix(uint32 i, EigenMatrix& Klocal);
		void ComputeLocalVolumetricStiffnessMatrixBB2(uint32 i, EigenMatrix& Klocal);

		void ComputeLocalJacobianBB2(uint32 elem, Vector3R v[10]) const; // not used

		// corotational helpers
		void ComputeCorotationalElasticForces(EigenVector& fout) const;
		void ComputeRotationMatrices();
		void AssembleDeviatoricStiffnessMatrixCR();
		void ComputeErrorCorotational(EigenVector& b, const SparseMatrix* J = nullptr) const;

		// explicit corotational
		void StepUnconstrainedCorotational(real h);
		void SolveConstraintsCorotationalSchur(real h);
		void SolveConstraintsCorotationalKKT(real h);

		// tractions
		bool NodeIsOnBoundary(int gidx) const { return mFixedPressureRows.find(gidx) != mFixedPressureRows.end(); }
		EigenVector ComputeTotalForce(bool corotational = false);

	protected:
		void AssembleJacobianMatrix(SparseMatrix& J, bool corot, bool update = false);
		void AssembleDeviatoricStiffnessMatrix();
		void ComputeLocalJacobian(uint32 e, Vector3R v[], bool update = false) const;

	protected:
		SparseMatrix mDeviatoricStiffnessMatrix; // deviatoric stiffness matrix
		EigenMatrix mBCStiffnessMatrix; // deviatoric stiffness matrix block for BCs
		EigenMatrix mVolumetricStiffnessMatrix; // volumetric stiffness matrix
		SparseMatrix mGeometricStiffnessMatrix;

		SparseMatrix mVolComplianceMatrix;
		SparseMatrix mVolJacobianMatrix;
		SparseMatrix mSparseSysMatrix;

		EigenMatrix mInverseMassMatrix;
		EigenMatrix mSystemMatrix;
		EigenMatrix mInverseImplicitMatrix;

		LinearSolver mSolver;
		
		uint32 mPressureOrder;
		uint32 mNumPressureNodes;
		std::unique_ptr<LinearTetrahedralMesh> mPressureMesh; // the pressure tet mesh

		EigenVector mElasticForce;
		std::vector<Matrix3R> mRotationMatrices;
		std::vector<EigenMatrix> mLocalStiffnessMatrices;

		Config mConfig;

		uint32 mNumFixedP; // number of prescribed pressure rows
		EigenMatrix mFixedJacobian; // Jacobian of prescribed pressure rows
		EigenMatrix mBCJacobian; // Jacobian of prescribed displacements
		std::vector<uint32> mPressureMap; // map from current (shuffled) pressure nodes to original ones
		std::vector<uint32> mInvPressureMap; // map from original pressure nodes to shuffled ones
		std::set<uint32> mFixedPressureRows; // nodes used for Neumann BCs

		bool mLogConstraint; // for nonlinear extension

		friend class FemTester;
	};

} // namespace FEM_SYSTEM

#endif // !FEM_PHYSICS_LINEAR_INCOMP_H
