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
#include "TetrahedralMesh.h"

namespace FEM_SYSTEM
{
	const float unit = 100.f;

	class FemPhysicsBase
	{
	public:
		FemPhysicsBase(const FemConfig& config);
		virtual ~FemPhysicsBase() {}

		virtual void Step(real dt) = 0;
		virtual Vector3R GetDeformedPosition(uint32 idx) const = 0;
		virtual uint32 GetNumNodes() const = 0;
		virtual void SolveEquilibrium(float) { }
		virtual void SetExternalForces(const std::vector<Vector3R>& forces) { }
		virtual void SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, const std::vector<uint32>& elemList, real pressure) { }
		virtual bool IsNodeFixed(uint32 i) const = 0;
		virtual TetrahedralMesh<uint32>* GetTetMesh() = 0;

		real GetLameFirstParam() const 
		{
			return real(mYoungsModulus * mPoissonRatio / (1.0 + mPoissonRatio) / (1.0 - 2.0 * mPoissonRatio));
		}

		real GetShearModulus() const // or Lame second param
		{
			return real(0.5 * mYoungsModulus / (1.0 + mPoissonRatio));
		}

		real GetBulkModulus() const
		{
			return real(mYoungsModulus / (1.0 - 2.0 * mPoissonRatio) / 3.0);
		}

	protected:
		MethodType mMethodType;
		real mYoungsModulus, mPoissonRatio;
		Matrix3R mNormalElasticityMatrix;
		Vector3R mGravity;
		int mNumSteps = 1;
		real mDensity; // body density
		SimulationType mSimType;

		// for applying external loads
		real mForceFraction = 0; // force application fraction
		real mForceStep = 0.1f;
	};

	inline FemPhysicsBase::FemPhysicsBase(const FemConfig& config)
		: mMethodType(config.mType)
		, mYoungsModulus(config.mYoungsModulus / unit)
		, mPoissonRatio(config.mPoissonRatio)
		, mGravity(0, config.mGravity * unit, 0)
		, mNumSteps(config.mNumSubsteps)
		, mDensity(config.mDensity / (unit * unit * unit))
		, mSimType(config.mSimType)
		, mForceStep(config.mForceApplicationStep)
	{
		// elastic params
		real omn = 1.f - mPoissonRatio;
		real om2n = 1.f - 2 * mPoissonRatio;
		real s = mYoungsModulus / (1.f + mPoissonRatio); // this is actually 2 * mu, i.e. twice the shear modulus
		real f = s / om2n;
		// the constitutive relation matrix for normal stresses/strains [Mueller]
		mNormalElasticityMatrix = f * Matrix3R(omn, mPoissonRatio, mPoissonRatio,
			mPoissonRatio, omn, mPoissonRatio,
			mPoissonRatio, mPoissonRatio, omn);
		// the constitutive relation matrix for shear components is just diag(2 * mu) [Mueller]
	}

} // namespace FEM_SYSTEM

#endif // FEM_PHYSICS_BASE_H
