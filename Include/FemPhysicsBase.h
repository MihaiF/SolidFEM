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

#ifndef FEM_PHYSICS_BASE_H
#define FEM_PHYSICS_BASE_H
#include "FemDataStructures.h"
#include "FemCollisions.h"
#include <Engine/Utils.h>

namespace FEM_SYSTEM
{
	class FemPhysicsBase
	{
	public:
		enum Axes
		{
			AXIS_X = 1,
			AXIS_Y = 2,
			AXIS_Z = 4,
		};

	public:
		FemPhysicsBase(const FemConfig& config);
		virtual ~FemPhysicsBase() {}

		virtual void Step(real dt) = 0;
		
		void HandleCollisions(real dt);

		// interface getters/setters
		virtual Vector3R GetDeformedPosition(uint32 idx) const = 0;
		virtual void SetDeformedPosition(uint32 idx, Vector3R val) { ASSERT(false) };
		virtual Vector3R GetInitialPosition(uint32 idx) const { ASSERT(false); return sZero; }
		virtual Vector3R GetVelocity(uint32 idx) const { ASSERT(false); return sZero; }
		virtual uint32 GetNumNodes() const = 0;
		virtual uint32 GetNumFreeNodes() const { ASSERT(false); return 0; }
		virtual uint32 GetNumLocalNodes() const { ASSERT(false); return 4; }
		virtual uint32 GetNumElements() const { ASSERT(false); return 0; }
		virtual void SetExternalForces(const std::vector<Vector3R>& forces) { ASSERT(false); }
		virtual void SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, real pressure) { ASSERT(false); }
		virtual bool IsNodeFixed(uint32 i) const = 0;
		virtual Matrix3R GetBarycentricJacobianMatrix(uint32 e) const { ASSERT(false); return Matrix3R(); }
		virtual real GetElementInitialVolume(uint32 e) const { ASSERT(false); return 0; }
		virtual uint32 GetGlobalIndex(uint32 e, uint32 l) const { ASSERT(false); return 0; }
		virtual uint32 GetGlobalOriginalIndex(uint32 e, uint32 l) const { ASSERT(false); return 0; }
		virtual real GetTotalVolume() const = 0;
		
		// interface methods
		virtual real ComputeEnergy(int level) const { ASSERT(false); return 0; }
		virtual void UpdatePositions(std::vector<Node>& nodes) { ASSERT(false); }
		virtual void ComputeDeformationGradient(uint32 e, Matrix3R& F) const { ASSERT(false); }

		// base class methods
		MaterialModelType GetMaterial() const { return mMaterial; }
		real GetYoungsModulus() const { return mYoungsModulus; }
		
		real GetLameFirstParam() const
		{
			return mLameLambda;
		}

		real GetShearModulus() const // or Lame second param
		{
			return real(0.5 * mYoungsModulus / (1.0 + mPoissonRatio));
		}

		real GetBulkModulus() const
		{
			return 1 / mInvBulkModulus;
		}

		real GetInverseBulkModulus()
		{
			return mInvBulkModulus;
		}

		real GetForceFraction() const { return mForceFraction; }

		void AddDirichletBC(uint32 i, uint32 axes);

		void GetNodeVolumetricStrains(std::vector<real>& volStrains) const;

		bool CheckForInversion(int verbose = 0);

		FemCollision* GetCollision() { return mFemCollision; }

	protected:
		void AssembleDynamicContributions();

	protected:
		MethodType mMethodType;
		real mYoungsModulus, mPoissonRatio;
		real mLameLambda;
		real mInvBulkModulus; // only used for mixed FEM
		Matrix3R mNormalElasticityMatrix;
		Vector3R mGravity;
		int mNumSteps = 1;
		real mDensity; // body density
		SimulationType mSimType;
		MaterialModelType mMaterial = MMT_LINEAR;
		uint32 mOuterIterations = 10;
		uint32 mInnerIterations = 100;
		real mContactStiffness;
		real mDirichletStiffness;
		
		// collisions
		bool mHasCollisions = false;
		FemCollision* mFemCollision;

		// Dirichlet BCs
		uint32 mNumBCs; // Dirichlet BC nodes
		std::vector<uint32> mReshuffleMap; // mapping from original nodes to shuffled ones

		// for applying external loads
		real mForceFraction = 0; // force application fraction
		real mForceStep = 0.1f;

		// variable Dirichlet BCs
		std::vector<uint32> mDirichletIndices;
		std::vector<uint32> mDirichletAxes;

		// Neumann/pressure BCs
		bool mUseImplicitPressureForces = false; // compute a stiffness matrix due to the variation of the normal
		std::vector<uint32> mTractionSurface; // triangle list of Neumann boundary
		real mAppliedPressure;
		EigenMatrix mTractionStiffnessMatrix;

		static Vector3R sZero;

		real mAbsNewtonResidualThreshold;

		bool mVerbose = false;

		SparseMatrix mContactJacobian;
		SparseMatrix mDirichletJacobian;

		real mTotalInitialVol;
	};

	inline void FemPhysicsBase::AddDirichletBC(uint32 i, uint32 axes)
	{
		uint32 idx = mReshuffleMap[i];
		mDirichletIndices.push_back(idx);
		mDirichletAxes.push_back(axes);
	}
} // namespace FEM_SYSTEM

#endif // FEM_PHYSICS_BASE_H
