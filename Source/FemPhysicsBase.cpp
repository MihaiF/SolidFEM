#include "FemPhysicsBase.h"
#include <Engine/Profiler.h>
#include "ElasticEnergy.h"

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
		mShearModulus = 0.5 * mYoungsModulus / (1.0 + mPoissonRatio);

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

	void FemPhysicsBase::AddCable(const std::vector<SpringNode>& cable, const Vector3Array& pos, real restLength, real stiffness, real damping, real actuation)
	{
		Cable femCable;
		femCable.mCableNodes = cable;
		femCable.mActuation = actuation;
		mCableRestLength = restLength * 0.01; // convert to meters
		mCableStiffness = stiffness;
		mCableDamping = damping;
		femCable.mCablePositions.resize(cable.size());
		for (uint32 i = 0; i < cable.size(); i++)
		{
			femCable.mCablePositions[i] = pos[i];
		}
		mCables.push_back(femCable);
	}

	void FemPhysicsBase::GetForceParamGrads(EigenVector& gradMu, EigenVector& gradLambda)
	{
		std::vector<Vector3R> fmu;
		std::vector<Vector3R> flambda;
		ElasticEnergy::ComputeForceParamGrads(this, fmu, flambda);

		gradMu = GetEigenVector(fmu, mNumBCs);
		gradLambda = GetEigenVector(flambda, mNumBCs);
	}

	void FemPhysicsBase::ComputeDeformationGradients()
	{
		int numElems = (int)GetNumElements();
		mDefGrads.resize(numElems);
		mDefGradsInv.resize(numElems);
		mDefGradsDet.resize(numElems);
		for (int e = 0; e < numElems; e++)
		{
			ComputeDeformationGradient(e, mDefGrads[e]);
			mDefGradsInv[e] = mDefGrads[e].GetInverse();
			mDefGradsDet[e] = mDefGrads[e].Determinant();
		}
	}

	real FemPhysicsBase::ComputeElasticEnergy()
	{
		real totalEnergy = 0;
		ComputeDeformationGradients();
		for (int e = 0; e < (int)GetNumElements(); e++)
		{
			totalEnergy += ElasticEnergy::ComputeElementEnergy(this, e);
		}
		return totalEnergy;
	}

	void FemPhysicsBase::ComputeSpringForces(Cable& cable, std::vector<Vector3R>& femForces)
	{
		if (cable.mCableNodes.empty() || cable.mActuation == 0)
			return;
		uint32 numNodes = (uint32)cable.mCableNodes.size();
		uint32 numSprings = numNodes - 1;
		// compute interpolated spring nodes
		cable.mCablePositions.resize(numNodes);
		Vector3Array vels(numNodes);
		for (uint32 i = 0; i < numNodes; i++)
		{
			int elem = cable.mCableNodes[i].elem;
			if (elem < 0)
			{
				cable.mCablePositions[i] = GetDeformedPosition(-elem);
				continue;
			}
			uint32 i0 = GetGlobalIndex(elem, 0);
			uint32 i1 = GetGlobalIndex(elem, 1);
			uint32 i2 = GetGlobalIndex(elem, 2);
			uint32 i3 = GetGlobalIndex(elem, 3);
			const Vector3R& x0 = GetDeformedPosition(i0, false);
			const Vector3R& x1 = GetDeformedPosition(i1, false);
			const Vector3R& x2 = GetDeformedPosition(i2, false);
			const Vector3R& x3 = GetDeformedPosition(i3, false);
			real w0 = cable.mCableNodes[i].bary.x;
			real w1 = cable.mCableNodes[i].bary.y;
			real w2 = cable.mCableNodes[i].bary.z;
			real w3 = 1 - w0 - w1 - w2;
			cable.mCablePositions[i] = w0 * x0 + w1 * x1 + w2 * x2 + w3 * x3;
			vels[i] = w0 * GetVelocity(i0) + w1 * GetVelocity(i1) + w2 * GetVelocity(i2) + w3 * GetVelocity(i3);
		}
		// 1. Compute spring errors and forces
		std::vector<Vector3R> springForces(numSprings);
		const real eps = mCableRestLength * 0.01;
		for (uint32 i = 0; i < numSprings; i++)
		{
			Vector3R y = cable.mCablePositions[i + 1] - cable.mCablePositions[i];
			real len = y.Length();
			Vector3R dir = (1.0 / len) * y;
			real err = len - mCableRestLength * cable.mActuation;
			real tension = -2 * mCableStiffness * err;
			// Bern tension model (unilateral spring + twice differentiable)
			if (err < -eps)
				tension = 0;
			else if (err >= -eps && err <= eps)
				tension = -mCableStiffness * (0.5 * err * err / eps + err + 0.5 * eps);

			// damping force (along spring)
			Vector3R z = vels[i + 1] - vels[i];
			real dv = z.Dot(dir);
			real damping = -mCableDamping * dv;

			springForces[i] = (tension + damping) * dir;
		}
		// 2. Compute spring node forces
		std::vector<Vector3R> nodeForces(numNodes);
		for (uint32 i = 0; i < numNodes; i++)
		{
			if (i < numSprings)
				nodeForces[i] -= springForces[i];
			if (i > 0)
				nodeForces[i] += springForces[i - 1];
		}
		// 3. Back-interpolate the FEM node forces
		for (uint32 i = 0; i < numNodes; i++)
		{
			int elem = cable.mCableNodes[i].elem;
			if (elem < 0)
			{
				femForces[-elem] = nodeForces[i];
				continue;
			}
			uint32 i0 = GetGlobalIndex(elem, 0);
			uint32 i1 = GetGlobalIndex(elem, 1);
			uint32 i2 = GetGlobalIndex(elem, 2);
			uint32 i3 = GetGlobalIndex(elem, 3);
			real w0 = cable.mCableNodes[i].bary.x;
			real w1 = cable.mCableNodes[i].bary.y;
			real w2 = cable.mCableNodes[i].bary.z;
			real w3 = 1 - w0 - w1 - w2;
			const Vector3R& f = mForceFraction * nodeForces[i];
			femForces[i0] += w0 * f;
			femForces[i1] += w1 * f;
			femForces[i2] += w2 * f;
			femForces[i3] += w3 * f;
		}
	}

	void FemPhysicsBase::ComputeSpringForces(std::vector<Vector3R>& femForces)
	{
		for (Cable& cable : mCables)
		{
			FemPhysicsBase::ComputeSpringForces(cable, femForces);
		}
	}

	template <class MATRIX>
	void FemPhysicsBase::ComputeSpringStiffnessMatrix(MATRIX& springStiffnessMatrix)
	{
		// finite difference implementation for now
		uint32 numNodes = GetNumFreeNodes();
		uint32 numDofs = numNodes * 3;
		springStiffnessMatrix.resize(numDofs, numDofs);
		springStiffnessMatrix.setZero();

		// current spring forces
		Vector3Array buffer(GetNumNodes());
		std::fill(buffer.begin(), buffer.end(), Vector3R());
		ComputeSpringForces(buffer);
		EigenVector fs0 = GetEigenVector(buffer, mNumBCs);

		real eps = 1e-5;
		for (uint32 j = 0; j < numDofs; j++)
		{
			Vector3R pos = GetDeformedPosition(j / 3 + mNumBCs, false);
			Vector3R inc;
			inc[j % 3] = eps;
			SetDeformedPosition(j / 3 + mNumBCs, pos + inc, false);
			std::fill(buffer.begin(), buffer.end(), Vector3R());
			ComputeSpringForces(buffer);
			EigenVector fs = GetEigenVector(buffer, mNumBCs);
			auto dfs = (1.0 / eps) * (fs - fs0);
			for (uint32 i = 0; i < numDofs; i++)
			{
				if (abs(dfs(i)) > eps)
					springStiffnessMatrix.coeffRef(i, j) = dfs(i);
			}
			SetDeformedPosition(j / 3 + mNumBCs, pos, false);
		}
	}


	template void FemPhysicsBase::ComputeSpringStiffnessMatrix(SparseMatrix& springStiffnessMatrix);
	template void FemPhysicsBase::ComputeSpringStiffnessMatrix(EigenMatrix& springStiffnessMatrix);
}