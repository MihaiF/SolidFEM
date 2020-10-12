#include "FemPhysicsBase.h"
#include <Engine/Profiler.h>

namespace FEM_SYSTEM
{
	const Vector3R sphPos(0.2f, -0.1f, 0);
	const real sphRad = 0.1f;

	Vector3R FemPhysicsBase::sZero;

	FemPhysicsBase::FemPhysicsBase(const FemConfig& config)
		: mMethodType(config.mType)
		, mYoungsModulus(config.mYoungsModulus)
		, mPoissonRatio(config.mPoissonRatio)
		, mNumSteps(config.mNumSubsteps)
		, mDensity(config.mDensity)
		, mSimType(config.mSimType)
		, mForceStep(config.mForceApplicationStep)
		, mMaterial(config.mMaterial)
		, mOuterIterations(config.mOuterIterations)
		, mInnerIterations(config.mInnerIterations)
		, mHasCollisions(config.mHasCollisions)
		, mContactStiffness(config.mContactStiffness)
		, mDirichletStiffness(1e4)
		, mAbsNewtonResidualThreshold(config.mAbsNewtonRsidualThreshold)
		, mVerbose(config.mVerbose)
		, mFemCollision(nullptr)
	{
		if (config.mZUpAxis)
			mGravity.Set(0, 0, config.mGravity);
		else
			mGravity.Set(0, config.mGravity, 0);

		if (config.mInvBulkModulus < 0)
			mInvBulkModulus = real((1.0 - 2.0 * mPoissonRatio) * 3.0 / mYoungsModulus); // default value
		else
			mInvBulkModulus = config.mInvBulkModulus;

		mLameLambda = real(mYoungsModulus * mPoissonRatio / (1.0 + mPoissonRatio) / (1.0 - 2.0 * mPoissonRatio));

		Printf("bulk modulus: %g\n", GetBulkModulus());
		Printf("shear modulus: %g\n", GetShearModulus());
		Printf("lambda: %g\n", GetLameFirstParam());

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

	void FemPhysicsBase::HandleCollisions(real dt)
	{
		PROFILE_SCOPE("Collision detection");
		if (!mHasCollisions)
			return;
		if (mFemCollision)
			mFemCollision->HandleCollisions(this);
	}

	void FemPhysicsBase::GetNodeVolumetricStrains(std::vector<real>& volStrains) const
	{
		std::vector<uint16> valences(GetNumNodes());
		for (uint32 e = 0; e < GetNumElements(); e++)
		{
			//auto E = ElasticEnergy::ComputeElementStrain(this, e);
			//real ev = E.Trace() / 3.f;
			Matrix3R F;
			ComputeDeformationGradient(e, F);
			real ev = F.Determinant() - 1; // J - 1
			for (uint32 i = 0; i < 4; i++)
			{
				uint32 globalIdx = GetGlobalIndex(e, i);
				valences[globalIdx]++;
				volStrains[globalIdx] += ev; // TODO: make sure with start with zero
			}
		}

		// average the element strains at nodes
		for (uint32 i = 0; i < GetNumNodes(); i++)
		{
			volStrains[i] *= 1.f / valences[i];
		}
	}

	bool FemPhysicsBase::CheckForInversion(int verbose)
	{
		bool ret = false;
		//if (verbose)
			//Printf("Check for inversion\n");
		for (uint32 i = 0; i < GetNumElements(); i++)
		{
			Matrix3R F;
			ComputeDeformationGradient(i, F);
			real det = F.Determinant();
			if (verbose > 1)
				Printf("%f\n", det);
			if (det <= 0.f)
			{
				ret = true;
				if (verbose > 1)
					Printf("Inverted element %d\n", i);
			}
		}

		if (ret && verbose)
			Printf("Inversion detected\n");

		return ret;
	}

	void FemPhysicsBase::AssembleDynamicContributions()
	{
		if (!mDirichletIndices.empty())
		{
			// count the constraints
			int count = 0;
			for (size_t i = 0; i < mDirichletIndices.size(); i++)
			{
				if (mDirichletAxes[i] & AXIS_X) count++;
				if (mDirichletAxes[i] & AXIS_Y) count++;
				if (mDirichletAxes[i] & AXIS_Z) count++;
			}

			mDirichletJacobian.resize(count, GetNumFreeNodes() * NUM_POS_COMPONENTS);
			mDirichletJacobian.setZero();
			count = 0;
			for (uint32 i = 0; i < mDirichletIndices.size(); i++)
			{
				uint32 idx = mDirichletIndices[i] - mNumBCs; // affecting node idx
				if (mDirichletAxes[i] & AXIS_X)
				{
					mDirichletJacobian.coeffRef(count, idx * 3) = 1;
					count++;
				}
				if (mDirichletAxes[i] & AXIS_Y)
				{
					mDirichletJacobian.coeffRef(count, idx * 3 + 1) = 1;
					count++;
				}
				if (mDirichletAxes[i] & AXIS_Z)
				{
					mDirichletJacobian.coeffRef(count, idx * 3 + 2) = 1;
					count++;
				}
			}
		}
	}


}