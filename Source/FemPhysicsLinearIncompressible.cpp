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

#include "FemPhysicsLinearIncompressible.h"
#include "MeshFactory.h"
#include "PolarDecomposition.h"
#include <Engine/Profiler.h>
#include <iostream>

#pragma warning( disable : 4244) // for double-float conversions
#pragma warning( disable : 4267) // for size_t-uint conversions

namespace FEM_SYSTEM
{
	FemPhysicsLinearIncompressible::FemPhysicsLinearIncompressible(std::vector<Tetrahedron>& tetrahedra,
		std::vector<Node>& nodes, 
		const FemConfig& config)
		: FemPhysicsLinear(tetrahedra, nodes, config)
		, mPressureOrder(1)
		, mNumFixedP(0)
	{
		if (config.mCustomConfig != nullptr)
		{
			mConfig = *((Config*)config.mCustomConfig);
			mPressureOrder = mConfig.mPressureOrder;
		}

		mNumPressureNodes = mPressureOrder == 0 ? GetNumElements() : nodes.size();

		// build the pressure tet mesh
		if (mPressureOrder > 0)
		{
			// use the linear mesh stored as an array of Tetrahedron in the base class
			StridedVector<uint16> stridedVec(&(tetrahedra[0].i[0]), tetrahedra.size(), sizeof(Tetrahedron));
			mPressureMesh.reset(MeshFactory::MakeMesh<TetrahedralMesh<uint32>>(mPressureOrder, nodes.size(), stridedVec, tetrahedra.size()));
		}

		AssembleDeviatoricStiffnessMatrix();

		// init pressure nodes maps to default (not shuffled)
		mPressureMap.resize(mNumPressureNodes);
		mInvPressureMap.resize(mNumPressureNodes);
		for (uint32 i = 0; i < mNumPressureNodes; i++)
		{
			mPressureMap[i] = i;
			mInvPressureMap[i] = i;
		}

		AssembleComplianceMatrix();
		AssembleJacobianMatrix(mVolJacobianMatrix, false);

		if (!mConfig.mUseKKTMatrix)
		{
			EigenMatrix denseM(mMassMatrix);
			mInverseMassMatrix = denseM.inverse();

			real h = 0.016; // FIXME
			EigenMatrix denseMK = EigenMatrix(mMassMatrix) + h * h * mDeviatoricStiffnessMatrix;
			mInverseImplicitMatrix = denseMK.inverse();
		}
	}

	// main stepping function
	void FemPhysicsLinearIncompressible::Step(real dt)
	{
		const int numSteps = 1; // small sub-steps for stability and less dissipation
		real h = dt / numSteps;
		for (int i = 0; i < numSteps; i++)
		{
			if (mSimType != ST_IMPLICIT)
			{
				if (!mConfig.mUseCorotational)
				{
					StepUnconstrained(h);
					if (!mConfig.mUseKKTMatrix)
						SolveConstraintsSchur(h);
					else
						SolveConstraintsKKT(h);
				}
				else
				{
					PROFILE_SCOPE("ICR explicit");
					ComputeRotationMatrices();
					StepUnconstrainedCorotational(h);
					if (!mConfig.mUseKKTMatrix)
						SolveConstraintsCorotationalSchur(h);
					else
						SolveConstraintsCorotationalKKT(h);
				}

			}
			else
			{
				if (mConfig.mUseCorotational)
				{
					ComputeRotationMatrices();
					if (mConfig.mUseKKTMatrix)
						StepImplicitCorotationalKKT(h);
					else
						StepImplicitCorotationalSchur(h);
				}
				else
				{
					if (mConfig.mUseKKTMatrix)
					{
						PROFILE_SCOPE("II-KKT");
						StepImplicitKKT(dt);
					}
					else
					{
						StepUnconstrained(h);
						SolveConstraintsSchur(h, true);
					}
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::SolveEquilibrium(float t)
	{
		EigenVector sol;
		if (mConfig.mUseStaticCondensation)
		{
			EigenMatrix C(mVolComplianceMatrix);
			EigenMatrix K = (mDeviatoricStiffnessMatrix + mVolJacobianMatrix.transpose() * C.inverse() * mVolJacobianMatrix); // static condensation
			auto decomp = K.fullPivLu();
			sol = decomp.solve(GetEigenVector(mBodyForces));
		}
		else
		{
			// for some reason broken for quadratic
			uint32 numNodes = GetNumFreeNodes();
			uint32 numPNodes = GetNumPressureNodes();
			uint32 numDofs = numNodes * 3;
			uint32 size = numDofs + numPNodes;

			mSystemMatrix.resize(size, size);
			mSystemMatrix.setZero();
			mSystemMatrix.block(0, 0, numDofs, numDofs) = mDeviatoricStiffnessMatrix;
			mSystemMatrix.block(0, numDofs, numDofs, numPNodes) = mVolJacobianMatrix.transpose();
			mSystemMatrix.block(numDofs, 0, numPNodes, numDofs) = mVolJacobianMatrix;
			if (!mConfig.mFullIncompressilble)
				mSystemMatrix.block(numDofs, numDofs, numPNodes, numPNodes) = -mVolComplianceMatrix;

			EigenVector rhs(size);
			rhs.setZero();
			rhs.head(numDofs) = GetEigenVector(mBodyForces);

			mSolver.Init(mSystemMatrix, LST_FULL_LU);
			EigenVector fullSol = mSolver.Solve(rhs);
			sol = fullSol.head(numDofs);
		}

		Vector3Array u = GetStdVector(sol);
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
			mDeformedPositions[i] += u[i];
	}

	void FemPhysicsLinearIncompressible::SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, const std::vector<uint32>& elemList, real pressure)
	{
		mAppliedPressure = pressure;
		mFixedPressureRows.clear();
		mTractionSurface.resize(triangleList.size());

		for (uint32 i = 0; i < triangleList.size(); i++)
		{
			int origIdx = triangleList[i]; // index in the original linear mesh
			int idx = mReshuffleMap[origIdx] - mNumBCs; // index in the reshuffled displacement nodes array
			ASSERT(idx >= 0);
			if (mPressureOrder == 1)
				mFixedPressureRows.insert(origIdx); // for pressure BCs (pressure per node)
			mTractionSurface[i] = idx;
		}

		if (mPressureOrder == 0)
		{
			for (uint32 i = 0; i < elemList.size(); i++)
				mFixedPressureRows.insert(elemList[i]);
		}
	}

	EigenVector FemPhysicsLinearIncompressible::ComputeTotalForce(bool corotational)
	{
		if (corotational)
			ComputeCorotationalElasticForces();
		else
		{
			uint32 numNodes = GetNumFreeNodes();

			// compute displacement vector
			Vector3Array u(numNodes);
			for (size_t i = 0; i < numNodes; i++)
			{
				u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
			}

			// compute explicit forces
			mElasticForce = mDeviatoricStiffnessMatrix * GetEigenVector(u);
		}

		EigenVector f = GetEigenVector(mBodyForces) - mElasticForce;
		{
			ComputeTractionForces();
			if (mTractionForces.size() == mBodyForces.size())
			{
				f += GetEigenVector(mTractionForces);
			}
		}
		if (mNumFixedP > 0)
		{
			EigenVector fixedP(mNumFixedP);
			fixedP.setZero();
			for (uint32 i = 0; i < mNumFixedP; i++)
				fixedP[i] = mAppliedPressure;
			EigenVector fp = mFixedJacobian.transpose() * fixedP;
			f -= fp;
		}

		if (mExternalForces.size() == mBodyForces.size())
		{
			f += GetEigenVector(mExternalForces);
		}

		return f;
	}

	void FemPhysicsLinearIncompressible::StepUnconstrained(real h)
	{
		EigenVector f = ComputeTotalForce();
		EigenVector sol = mInverseMassMatrix * f;
		Vector3Array aext = GetStdVector(sol);

		// candidate positions - these are only needed for the strain (in the RHS)
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
		{
			mVelocities[i] += h * aext[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsLinearIncompressible::SolveConstraintsSchur(real h, bool implicit)
	{
		uint32 numNodes = GetNumFreeNodes();

		// compute displacement vector
		Vector3Array u(numNodes);
		for (size_t i = 0; i < numNodes; i++)
		{
			u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
		}

		const auto& J = mVolJacobianMatrix;
		auto Jt = J.transpose();
		const auto& Minv = implicit ? mInverseImplicitMatrix : mInverseMassMatrix;
		EigenVector b = J * GetEigenVector(u);
		EigenVector a;
		if (implicit)
		{
			a = h * h * mDeviatoricStiffnessMatrix * GetEigenVector(mVelocities);
			b -= h * J * mInverseImplicitMatrix * a;
		}

		// form the linear system and solve it
		if (mSystemMatrix.rows() == 0)
		{
			mSystemMatrix = h * h * J * Minv * Jt;
			if (!mConfig.mFullIncompressilble)
				mSystemMatrix += mVolComplianceMatrix;
			mSolver.Init(mSystemMatrix, LST_LLT);
		}
		EigenVector lambda = mSolver.Solve(b);
		EigenVector f = h * Jt * lambda;
		if (implicit)
			f += a;

		// solve the linear system for the accelerations
		EigenVector sol = Minv * f;
		Vector3Array acc = GetStdVector(sol);

		// integrate node velocities and positions
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] -= acc[i];
			mDeformedPositions[i] -= h * acc[i];
		}
	}

	// experimental
	void SolveUzawa(const EigenMatrix& A, const EigenMatrix& J, const EigenMatrix& C, 
		const EigenVector& a, const EigenVector& b, EigenVector& dv)
	{
		uint32 numConstr = J.rows();
		EigenVector p(numConstr);
		p.setZero();
		auto decomp = A.lu();
		const real omega = 0.1f;
		for (uint32 iter = 0; iter < 100; iter++)
		{
			dv = decomp.solve(a - J.transpose() * p);
			p += omega * (J * dv - C * p - b);
		}
	}

	void FemPhysicsLinearIncompressible::SolveConstraintsKKT(real h, bool implicit)
	{
		PROFILE_SCOPE("Incomp explicit KKT");
		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumFreePressureNodes();

		Vector3Array u(numNodes);
		for (size_t i = 0; i < numNodes; i++)
		{
			u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
		}
		EigenVector disp = GetEigenVector(u);
		EigenVector b = -mVolJacobianMatrix * disp;

		uint32 numDofs = numNodes * 3;
		uint32 size = numDofs + numPNodes;
		EigenVector rhs(size);
		rhs.setZero();
		EigenVector a(numDofs);
		a.setZero();
		if (implicit)
		{
			auto K = mDeviatoricStiffnessMatrix;
			a = -h * h * K * GetEigenVector(mVelocities);
			rhs.head(numDofs) = a;
		}
		rhs.tail(numPNodes) = b;

		// form the linear system and solve it
		if (mSystemMatrix.rows() == 0)
		{
			mSystemMatrix.resize(size, size);
			mSystemMatrix.setZero();
			EigenMatrix M = mMassMatrix;
			if (implicit)
			{
				M += h * h * mDeviatoricStiffnessMatrix;
			}
			mSystemMatrix.block(0, 0, numDofs, numDofs) = M;
			mSystemMatrix.block(0, numDofs, numDofs, numPNodes) = h * mVolJacobianMatrix.transpose();
			mSystemMatrix.block(numDofs, 0, numPNodes, numDofs) = h * mVolJacobianMatrix;
			if (!mConfig.mFullIncompressilble)
				mSystemMatrix.block(numDofs, numDofs, numPNodes, numPNodes) = -mVolComplianceMatrix;
			mSolver.Init(mSystemMatrix, LST_LDLT);
		}
	
		EigenVector dv, p;
		EigenVector sol = mSolver.Solve(rhs);
		dv = sol.head(numDofs);
		p = sol.tail(numPNodes);

		Vector3Array deltas = GetStdVector(dv);
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] += deltas[i];
			mDeformedPositions[i] += h * deltas[i];
		}
	}

	void FemPhysicsLinearIncompressible::StepImplicitKKT(real h)
	{
		PROFILE_SCOPE("II-KKT");
		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumFreePressureNodes();

		Vector3Array u(numNodes);
		for (size_t i = 0; i < numNodes; i++)
		{
			u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
		}
		EigenVector uEigen = GetEigenVector(u);
		EigenVector b = -mVolJacobianMatrix * uEigen;

		EigenVector f = ComputeTotalForce();
		EigenVector a = mMassMatrix * GetEigenVector(mVelocities) + h * f;

		uint32 numDofs = numNodes * 3;
		uint32 size = numDofs + numPNodes;
		EigenVector rhs(size);
		rhs.head(numDofs) = a;
		rhs.tail(numPNodes) = b;

		// form the linear system and solve it
		if (mSystemMatrix.rows() == 0)
		{
			mSystemMatrix.resize(size, size);
			mSystemMatrix.setZero();
			mSystemMatrix.block(0, 0, numDofs, numDofs) = EigenMatrix(mMassMatrix) + h * h * mDeviatoricStiffnessMatrix;
			mSystemMatrix.block(0, numDofs, numDofs, numPNodes) = h * mVolJacobianMatrix.transpose();
			mSystemMatrix.block(numDofs, 0, numPNodes, numDofs) = h * mVolJacobianMatrix;
			if (!mConfig.mFullIncompressilble)
				mSystemMatrix.block(numDofs, numDofs, numPNodes, numPNodes) = -mVolComplianceMatrix;
			mSolver.Init(mSystemMatrix, LST_LDLT, true);
		}
		EigenVector sol = mSolver.Solve(rhs);

		EigenVector v = sol.head(numDofs);
		Vector3Array vels = GetStdVector(v);
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] = vels[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	// old version solving directly for final velocities (not for deltas, without unconstrained step)
	void FemPhysicsLinearIncompressible::StepImplicitSchur(real h)
	{
		PROFILE_SCOPE("II-Schur");
		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumPressureNodes();

		Vector3Array u(numNodes);
		for (size_t i = 0; i < numNodes; i++)
		{
			u[i] = mDeformedPositions[i] - mReferencePositions[i + mNumBCs];
		}
		EigenVector uEigen = GetEigenVector(u);
		EigenVector b = -mVolJacobianMatrix * uEigen;
		EigenVector a = mMassMatrix * GetEigenVector(mVelocities) - h * mDeviatoricStiffnessMatrix * uEigen + h * GetEigenVector(mBodyForces);

		// form the linear system and solve it
		if (mSystemMatrix.rows() == 0)
		{
			mSystemMatrix = h * h * mVolJacobianMatrix * mInverseImplicitMatrix * mVolJacobianMatrix.transpose() + EigenMatrix(mVolComplianceMatrix);
			mSolver.Init(mSystemMatrix, LST_LDLT);
		}
		EigenVector rhs = -b + h * mVolJacobianMatrix * mInverseImplicitMatrix * a;
		EigenVector p = mSolver.Solve(rhs);
		EigenVector v = mInverseImplicitMatrix * (a - h * mVolJacobianMatrix.transpose() * p);
		Vector3Array vels = GetStdVector(v);
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] = vels[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsLinearIncompressible::StepImplicitCorotationalKKT(real h)
	{
		PROFILE_SCOPE("IICR-KKT");

		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumFreePressureNodes();

		EigenMatrix J;
		AssembleJacobianMatrix(J, true);

		EigenVector b;
		ComputeErrorCorotational(b);

		EigenVector f = ComputeTotalForce(true);
		EigenVector a = mMassMatrix * GetEigenVector(mVelocities) + h * f;

		real flip = 1; // solvers that require symmetric matrices like LDLT or LLT are excluded for -1 (only LU left)

		uint32 numDofs = numNodes * 3;
		uint32 size = numDofs + numPNodes;
		EigenVector rhs(size);
		rhs.head(numDofs) = a;
		rhs.tail(numPNodes) = -b;

		// form the linear system and solve it
		AssembleDeviatoricStiffnessMatrixCR();
		EigenVector sol;
		{
			PROFILE_SCOPE("IICR build S");
			if (mSystemMatrix.rows() == 0)
			{
				mSystemMatrix.resize(size, size);
				mSystemMatrix.block(numDofs, numDofs, numPNodes, numPNodes) = -flip * mVolComplianceMatrix;
			}
			mSystemMatrix.block(0, 0, numDofs, numDofs) = EigenMatrix(mMassMatrix)
				+ h * h * mDeviatoricStiffnessMatrix.block(mNumBCs * 3, mNumBCs * 3, numDofs, numDofs);
			mSystemMatrix.block(0, numDofs, numDofs, numPNodes) = flip * h * J.transpose();
			mSystemMatrix.block(numDofs, 0, numPNodes, numDofs) = h * J;
		}
		{
			PROFILE_SCOPE("IICR solve");
			mSolver.Init(mSystemMatrix, LST_LDLT, true);
			sol = mSolver.Solve(rhs);
		}

		EigenVector v = sol.head(numDofs);
		Vector3Array vels = GetStdVector(v);
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] = vels[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsLinearIncompressible::StepImplicitCorotationalSchur(real h)
	{
		PROFILE_SCOPE("IICR-Schur");

		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumPressureNodes();

		EigenVector b;
		ComputeErrorCorotational(b);

		ComputeCorotationalElasticForces();
		EigenVector a = mMassMatrix * GetEigenVector(mVelocities) - h * mElasticForce + h * GetEigenVector(mBodyForces);

		AssembleDeviatoricStiffnessMatrixCR();
		uint32 numDofs = numNodes * 3;
		EigenMatrix Kcr = mDeviatoricStiffnessMatrix.block(mNumBCs * 3, mNumBCs * 3, numDofs, numDofs);

		// form the linear system and solve it
		BEGIN_PROFILE("Form lin sys");
		SparseMatrix J;
		AssembleJacobianMatrix(J, true);
		auto Jt = J.transpose();
		EigenMatrix M = mMassMatrix + h * h * Kcr;
		EigenMatrix Minv = M.inverse();
		EigenMatrix JMinv = J * Minv;
		EigenVector rhs = b + h * JMinv * a;
		EigenMatrix S = h * h * JMinv * Jt + EigenMatrix(mVolComplianceMatrix);
		END_PROFILE();

		BEGIN_PROFILE("Solve")
		auto solver = S.llt();
		
		EigenVector p = solver.solve(rhs);
		END_PROFILE();
		
		EigenVector v = Minv * (a - h * Jt * p);
		Vector3Array vels = GetStdVector(v);
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] = vels[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsLinearIncompressible::ComputeRotationMatrices()
	{
		PROFILE_SCOPE("Compute rot");
		mRotationMatrices.resize(GetNumElements());
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			// compute the deformation gradient and its polar decomposition
			Matrix3R F, U;
			ComputeDeformationGradient(i, F);
			ComputePolarDecomposition(F, mRotationMatrices[i], U);
		}
	}

	void FemPhysicsLinearIncompressible::ComputeCorotationalElasticForces()
	{
		size_t numNodes = GetNumFreeNodes(); // remove the first 4 entries which are fixed (hack)
		size_t numDofs = numNodes * 3;
		mElasticForce.resize(numDofs, 1);
		mElasticForce.setZero();

		// go through all elements (tetrahedra)
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			EigenMatrix Klocal = mLocalStiffnessMatrices[i];

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
					f += Kblock * (!mRotationMatrices[i] * GetDeformedPosition(mTetMesh->GetGlobalIndex(i, k)) - mReferencePositions[kGlobal]);
				}

				// rotate the force and add it to the global one
				f = mRotationMatrices[i] * f;
				for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
				{
					mElasticForce[jOffset + x] += f[x];
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::AssembleDeviatoricStiffnessMatrix()
	{
		size_t numNodes = GetNumFreeNodes(); // remove the first 4 entries which are fixed (hack)
		size_t numDofs = numNodes * 3;
		mDeviatoricStiffnessMatrix.resize(numDofs, numDofs);
		mDeviatoricStiffnessMatrix.setZero();
		mLocalStiffnessMatrices.resize(GetNumElements());
		// go through all linear elements (tetrahedra)
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			EigenMatrix Klocal;
			real G = GetShearModulus();
			ComputeLocalStiffnessMatrixBB(i, G, -2.0 / 3.0 * G, Klocal);
			mLocalStiffnessMatrices[i] = Klocal;
			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (size_t j = 0; j < GetNumLocalNodes(); j++)
			{
				size_t jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, j)];
				for (size_t k = 0; k < GetNumLocalNodes(); k++)
				{
					size_t kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, k)];
					if (jGlobal < mNumBCs || kGlobal < mNumBCs)
						continue;
					int jOffset = (jGlobal - mNumBCs) * NUM_POS_COMPONENTS;
					int kOffset = (kGlobal - mNumBCs) * NUM_POS_COMPONENTS;

					// add the the whole 3x3 block to the block matrix
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							mDeviatoricStiffnessMatrix.coeffRef(jOffset + x, kOffset + y) += Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
						}
					}
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::ComputeLocalDeviatoricStiffnessMatrix(uint32 i, EigenMatrix& Klocal)
	{
		// this code only works for linear elements
		ASSERT(GetNumLocalNodes() == 4);
		Klocal = EigenMatrix(12, 12);

		Matrix3R Bn[4], Bs[4];
		ComputeStrainJacobian(i, Bn, Bs);

		// project the normal matrices
		Matrix3R P(2. / 3., -1. / 3., -1. / 3.,
			-1. / 3., 2. / 3., -1. / 3.,
			-1. / 3., -1. / 3., 2. / 3.);
		for (int j = 0; j < 4; j++)
			Bn[j] = P * Bn[j];
		// TODO: reuse code

		// for linear FEM we can actually precompute the tangent stiffness matrix
		real shearModulus = 0.5f * mYoungsModulus / (1.f + mPoissonRatio);
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				Matrix3R block = mElementVolumes[i] * shearModulus * (2.0 * !Bn[j] * Bn[k] + Bs[j] * Bs[k]);
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

	void FemPhysicsLinearIncompressible::ComputeLocalVolumetricStiffnessMatrix(uint32 i, EigenMatrix& Klocal)
	{
		// this code only works for linear elements
		ASSERT(GetNumLocalNodes() == 4);
		Klocal = EigenMatrix(12, 12);

		// for linear FEM we can actually precompute the tangent stiffness matrix
		real bulkModulus = mYoungsModulus / (1. - 2. * mPoissonRatio) / 3.;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				Matrix3R L = Matrix3R::TensorProduct(mBarycentricJacobians[i].y[j], mBarycentricJacobians[i].y[k]);
				Matrix3R block = mElementVolumes[i] * bulkModulus * L;
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

	void FemPhysicsLinearIncompressible::AssembleComplianceMatrix()
	{
		//MEASURE_TIME("Assemble C");
		real bulkModulus = GetBulkModulus();
		real Cconst = 1.f / bulkModulus;

		// the local compliance matrix template
		const uint32 numLocalPressureNodes = GetNumLocalPressureNodes();
		const uint32 numLocalPressureDofs = numLocalPressureNodes * NUM_STRESS_COMPONENTS;
		EigenMatrix Clocal(numLocalPressureDofs, numLocalPressureDofs);
		Clocal.setZero();
		for (uint32 i = 0; i < numLocalPressureNodes; i++)
		{
			for (uint32 j = 0; j < numLocalPressureNodes; j++)
			{
				real factor = 1.f;
				if (mPressureOrder == 1)
				{
					factor = i == j ? 0.1 : 0.05;
				}
				Clocal(i, j) = factor * Cconst;
			}
		}

		// assemble the global compliance matrix 
		uint32 numPressureNodes = GetNumFreePressureNodes();
		mVolComplianceMatrix = SparseMatrix(numPressureNodes, numPressureNodes);
		mVolComplianceMatrix.setZero();
		for (uint32 e = 0; e < GetNumElements(); e++)
		{
			real vol = mElementVolumes[e];
			for (uint32 i = 0; i < numLocalPressureNodes; i++)
			{
				uint32 globalI = mInvPressureMap[GetPressureGlobalIndex(e, i)];
				if (globalI >= numPressureNodes)
					continue; // skip fixed pressure values
				for (uint32 j = 0; j < numLocalPressureNodes; j++)
				{
					uint32 globalJ = mInvPressureMap[GetPressureGlobalIndex(e, j)];
					if (globalJ >= numPressureNodes)
						continue; // skip fixed pressure values
					mVolComplianceMatrix.coeffRef(globalI, globalJ) += vol * Clocal(i, j);
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::ComputeLocalJacobian(uint32 e, Vector3R v[])
	{
		if (mOrder == 1)
		{
			for (int i = 0; i < 4; i++)
				v[i] = mBarycentricJacobians[e].y[i];
		}
		else
			ComputeLocalJacobianBB2(e, v);
	}

	template<class MATRIX>
	void FemPhysicsLinearIncompressible::AssembleJacobianMatrix(MATRIX& J, bool corot)
	{
		PROFILE_SCOPE("IICR assemble J");
		uint32 numNodes = GetNumFreeNodes();
		uint32 nDof = NUM_POS_COMPONENTS * numNodes; // position DOFs
		uint32 numLocalNodes = GetNumLocalNodes();
		uint32 nLocalDofs = numLocalNodes * NUM_POS_COMPONENTS;

		uint32 numPressureLocalNodes = GetNumLocalPressureNodes();
		uint32 numPNodes = GetNumFreePressureNodes();

		J.resize(numPNodes, nDof);
		J.setZero();
		mFixedJacobian = EigenMatrix(mNumFixedP, nDof);
		mFixedJacobian.setZero();
		real factor = mPressureOrder == 1 ? 0.25f : 1.f;
		int fixedCounter = 0;
		for (uint32 e = 0; e < GetNumElements(); e++)
		{
			real vol = mElementVolumes[e];
			Vector3R v[10];
			ComputeLocalJacobian(e, v);

			if (corot)
			{
				for (uint32 k = 0; k < numLocalNodes; k++)
					v[k] = mRotationMatrices[e] * v[k];
			}

			EigenMatrix Jlocal(1, nLocalDofs);
			Jlocal.setZero();
			{
				// build Jacobian matrix
				for (uint32 k = 0; k < numLocalNodes; k++)
				{
					for (int l = 0; l < NUM_POS_COMPONENTS; l++)
					{
						Jlocal(0, k * NUM_POS_COMPONENTS + l) = v[k][l];
					}
				}
			}

			// go through all local stress nodes
			for (uint32 j = 0; j < numPressureLocalNodes; j++)
			{
				uint32 globalJ = mInvPressureMap[GetPressureGlobalIndex(e, j)];
				// go through all local position nodes
				for (uint32 k = 0; k < GetNumLocalNodes(); k++)
				{
					uint32 globalK = mReshuffleMap[mTetMesh->GetGlobalIndex(e, k)];
					if (globalK < mNumBCs)
						continue;
					int baseCol = (globalK - mNumBCs) * NUM_POS_COMPONENTS;
					if (mOrder == 2 && mPressureOrder == 1)
					{
						MultiIndex C = mPressureMesh->mIJKL + j * 4;
						Vector3R v;
						for (int n = 0; n < 4; n++)
						{
							MultiIndex A = mTetMesh->mIJKL + k * 4;
							if (A.Decrement(n))
								v += 0.2f * ComputeMultiIndexSumFactor(1, C, A) * mBarycentricJacobians[e].y[n];
						}
						if (corot)
							v = mRotationMatrices[e] * v;
						for (int l = 0; l < NUM_POS_COMPONENTS; l++)
							J.coeffRef(globalJ, baseCol + l) += vol * v[l];
					}
					else
					{
						for (int l = 0; l < NUM_POS_COMPONENTS; l++)
						{
							if (globalJ >= numPNodes)
								mFixedJacobian.coeffRef(globalJ - numPNodes, baseCol + l) += factor * vol * Jlocal(0, k * NUM_POS_COMPONENTS + l);
							else
								J.coeffRef(globalJ, baseCol + l) += factor * vol * Jlocal(0, k * NUM_POS_COMPONENTS + l);
						}
					}
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::ComputeLocalDeviatoricStiffnessMatrixBB2(uint32 k, EigenMatrix& Klocal)
	{
		ASSERT(GetNumLocalNodes() == 10);
		Klocal = EigenMatrix(30, 30);

		real mu = 0.5f * mYoungsModulus / (1.f + mPoissonRatio); // shear modulus

		Matrix3R Bn[4], Bs[4];
		ComputeStrainJacobian(k, Bn, Bs);

		// project the normal matrices
		Matrix3R P(2.f / 3.f, -1.f / 3.f, -1.f / 3.f,
			-1.f / 3.f, 2.f / 3.f, -1.f / 3.f,
			-1.f / 3.f, -1.f / 3.f, 2.f / 3.f);
		for (int j = 0; j < 4; j++)
			Bn[j] = P * Bn[j];

		Matrix3R identity;
		for (uint32 i = 0; i < GetNumLocalNodes(); i++)
		{
			for (uint32 j = 0; j < GetNumLocalNodes(); j++)
			{
				int* multiIndexI = mTetMesh->mIJKL + i * 4;
				int* multiIndexJ = mTetMesh->mIJKL + j * 4;
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
						Kn = Kn + factor * !Bn[c] * Bn[d];
						Ks = Ks + factor * Bs[c] * Bs[d];
					}
				}
				Matrix3R block = 0.4 * mElementVolumes[k] * mu * (2.f * Kn + Ks);

				// copy the block into Klocal
				for (uint32 a = 0; a < NUM_POS_COMPONENTS; a++)
				{
					for (uint32 b = 0; b < NUM_POS_COMPONENTS; b++)
					{
						Klocal(i * NUM_POS_COMPONENTS + a, j * NUM_POS_COMPONENTS + b) = 4 * block(a, b);
					}
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::ComputeLocalVolumetricStiffnessMatrixBB2(uint32 k, EigenMatrix& Klocal)
	{
		ASSERT(GetNumLocalNodes() == 10);
		Klocal = EigenMatrix(30, 30);

		real bulkModulus = mYoungsModulus / (1.f - 2.f * mPoissonRatio) / 3.f;

		Matrix3R identity;
		const auto& y = mBarycentricJacobians[k].y;
		for (uint32 i = 0; i < GetNumLocalNodes(); i++)
		{
			for (uint32 j = 0; j < GetNumLocalNodes(); j++)
			{
				int* multiIndexI = mTetMesh->mIJKL + i * 4;
				int* multiIndexJ = mTetMesh->mIJKL + j * 4;
				Matrix3R Kn(0);
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
						Kn = Kn + factor * Matrix3R::TensorProduct(y[c], y[d]);
					}
				}
				Matrix3R block = 0.4f * mElementVolumes[k] * bulkModulus * Kn;

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

	void FemPhysicsLinearIncompressible::ComputeLocalJacobianBB2(uint32 elem, Vector3R v[10])
	{
		for (int A = 0; A < 10; A++)
		{
			v[A].SetZero();
			for (int n = 0; n < 4; n++)
			{
				MultiIndex B = mTetMesh->mIJKL + A * 4;
				if (B.Decrement(n))
					v[A] += mBarycentricJacobians[elem].y[n];
			}
		}
	}

	void FemPhysicsLinearIncompressible::AssembleDeviatoricStiffnessMatrixCR()
	{
		PROFILE_SCOPE("IICR assemble K");
		size_t numNodes = GetNumNodes();
		size_t numDofs = numNodes * 3;
		if (mDeviatoricStiffnessMatrix.rows() != numDofs)
			mDeviatoricStiffnessMatrix.resize(numDofs, numDofs);
		mDeviatoricStiffnessMatrix.setZero();
		// go through all tetrahedra
		for (size_t i = 0; i < GetNumElements(); i++)
		{
			EigenMatrix Klocal = mLocalStiffnessMatrices[i];

			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (size_t j = 0; j < GetNumLocalNodes(); j++)
			{
				size_t jGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, j)];
				for (size_t k = 0; k < GetNumLocalNodes(); k++)
				{
					size_t kGlobal = mReshuffleMap[mTetMesh->GetGlobalIndex(i, k)];
					int jOffset = jGlobal * NUM_POS_COMPONENTS;
					int kOffset = kGlobal * NUM_POS_COMPONENTS;

					// construct the 3x3 block
					Matrix3R Kblock;
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							Kblock(x, y) = Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
						}
					}

					const Matrix3R& R = mRotationMatrices[i];
					Matrix3R Krot = R * Kblock * !R;

					// add the the whole 3x3 block to the block matrix
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							mDeviatoricStiffnessMatrix.coeffRef(jOffset + x, kOffset + y) += Krot(x, y);
						}
					}
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::StepUnconstrainedCorotational(real h)
	{
		ComputeCorotationalElasticForces();
		EigenVector sol = mInverseMassMatrix * (GetEigenVector(mBodyForces) - mElasticForce);
		Vector3Array aext = GetStdVector(sol);

		// candidate positions - these are only needed for the strain (in the RHS)
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
		{
			mVelocities[i] += h * aext[i];
			mDeformedPositions[i] += h * mVelocities[i];
		}
	}

	void FemPhysicsLinearIncompressible::ComputeErrorCorotational(EigenVector& b, const EigenMatrix* J)
	{
		if (J != nullptr)
		{
			uint32 numNodes = GetNumFreeNodes();
			Vector3Array x(numNodes);
			Vector3Array X(numNodes);
			for (size_t i = 0; i < numNodes; i++)
			{
				x[i] = mDeformedPositions[i];
				X[i] = mReferencePositions[i + mNumBCs];
			}
			b = *J * GetEigenVector(x) - mVolJacobianMatrix * GetEigenVector(X);
		}
		else
		{
			uint32 numPressureLocalNodes = GetNumLocalPressureNodes();
			uint32 numPNodes = GetNumFreePressureNodes();
				
			b.resize(numPNodes);
			b.setZero();
			for (uint32 e = 0; e < GetNumElements(); e++)
			{
				Vector3R v[10];
				ComputeLocalJacobian(e, v);
				real factor = mPressureOrder == 1 ? 0.25f : 1.f;

				// compute displacements and their contribution
				real sum = 0;
				Vector3R disp[10];
				for (uint32 j = 0; j < GetNumLocalNodes(); j++)
				{
					uint32 globalIdx = mTetMesh->GetGlobalIndex(e, j);
					uint32 globalIdxMapped = mReshuffleMap[globalIdx];
					disp[j] = !mRotationMatrices[e] * GetDeformedPosition(globalIdx) - mReferencePositions[globalIdxMapped];
					sum += factor * mElementVolumes[e] * dot(v[j], disp[j]);
				}

				for (uint32 j = 0; j < numPressureLocalNodes; j++)
				{
					uint32 globalJ = mInvPressureMap[GetPressureGlobalIndex(e, j)];
					if (globalJ >= numPNodes)
						continue; // skip fixed pressure nodes
					if (mOrder == 2 && mPressureOrder == 1)
					{
						sum = 0;
						for (uint32 k = 0; k < GetNumLocalNodes(); k++)
						{
							MultiIndex C = mPressureMesh->mIJKL + j * 4;
							Vector3R v;
							for (int n = 0; n < 4; n++)
							{
								MultiIndex A = mTetMesh->mIJKL + k * 4;
								if (A.Decrement(n))
									v += 0.2f * ComputeMultiIndexSumFactor(1, C, A) * mBarycentricJacobians[e].y[n];
							}							
							sum += mElementVolumes[e] * dot(v, disp[k]);
						}
					}
					b[globalJ] += sum;
				}
			}
		}
	}

	void FemPhysicsLinearIncompressible::SolveConstraintsCorotationalSchur(real h)
	{
		PROFILE_SCOPE("ICR-Schur");
		uint32 numNodes = GetNumFreeNodes();

		// form the linear system and solve it
		EigenMatrix J;
		AssembleJacobianMatrix(J, true);
		auto Jt = J.transpose();

		EigenVector b;
		{
			PROFILE_SCOPE("ICR-Schur build sys");
			mSystemMatrix = h * h * J * mInverseMassMatrix * Jt + EigenMatrix(mVolComplianceMatrix);

			ComputeErrorCorotational(b);
		}

		EigenVector lambda;
		{
			PROFILE_SCOPE("ICR-Schur solve");
			mSolver.Init(mSystemMatrix, LST_LLT);
			lambda = mSolver.Solve(b);
		}
		EigenVector f = Jt * lambda;

		// solve the linear system for the accelerations
		EigenVector sol = mInverseMassMatrix * (f);
		Vector3Array a = GetStdVector(sol);

		// integrate node velocities and positions
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] -= h * a[i];
			mDeformedPositions[i] -= h * h * a[i];
		}
	}

	void FemPhysicsLinearIncompressible::SolveConstraintsCorotationalKKT(real h)
	{
		PROFILE_SCOPE("Solve CR KKT");
		uint32 numNodes = GetNumFreeNodes();
		uint32 numPNodes = GetNumPressureNodes();

		EigenMatrix J;
		AssembleJacobianMatrix(J, true);

		EigenVector b;
		ComputeErrorCorotational(b);

		uint32 numDofs = numNodes * 3;
		uint32 size = numDofs + numPNodes;
		EigenVector rhs(size);
		rhs.setZero();
		rhs.tail(numPNodes) = -b;

		// form the linear system and solve it
		if (mSystemMatrix.rows() == 0)
		{
			mSystemMatrix.resize(size, size);
			mSystemMatrix.setZero();
			mSystemMatrix.block(0, 0, numDofs, numDofs) = mMassMatrix;
			if (!mConfig.mFullIncompressilble)
				mSystemMatrix.block(numDofs, numDofs, numPNodes, numPNodes) = -mVolComplianceMatrix;
		}
		
		mSystemMatrix.block(0, numDofs, numDofs, numPNodes) = h * J.transpose();
		mSystemMatrix.block(numDofs, 0, numPNodes, numDofs) = h * J;

		mSolver.Init(mSystemMatrix, LST_LDLT, true);
		EigenVector sol = mSolver.Solve(rhs);

		EigenVector v = sol.head(numDofs);
		Vector3Array vels = GetStdVector(v);
		for (size_t i = 0; i < numNodes; i++)
		{
			mVelocities[i] += vels[i];
			mDeformedPositions[i] += h * vels[i];
		}
	}

} // namespace FEM_SYSTEM