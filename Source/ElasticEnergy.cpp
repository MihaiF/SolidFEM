/*
BSD 3-Clause License

Copyright (c) 2020, Mihai Francu
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

#include "ElasticEnergy.h"
#include "FemPhysicsBase.h"
#include "PolarDecomposition.h"
#include <Engine/Profiler.h>

namespace FEM_SYSTEM
{
	template <typename TE, typename TR>
	TR InnerProduct(const std::vector<TE>& a, const std::vector<TE>& b)
	{
		ASSERT(a.size() == b.size());
		TR sum = 0;
		for (size_t i = 0; i < a.size(); i++)
			sum += a[i] * b[i];
		return sum;
	}

	Matrix3R ElasticEnergy::ComputeElementStress(const FemPhysicsBase* femPhysics, uint32 e, int material)
	{
		Matrix3R id;
		const real mu = femPhysics->GetShearModulus();
		const real lambda = femPhysics->GetLameFirstParam();
		if (material == -1)
			material = femPhysics->GetMaterial();

		// compute nodal forces for this element only
		Matrix3R F;
		femPhysics->ComputeDeformationGradient(e, F);

		Matrix3R R, U;
		if (material == MMT_COROTATIONAL)
		{
			ComputePolarDecomposition(F, R, U);
			F = U; // this seems to be broken for the torus!
		}

		// compute strain and strain rate as tensors
		Matrix3R strain;
		if (material == MMT_LINEAR ||
			material == MMT_COROTATIONAL ||
			material == MMT_DISTORTIONAL_LINEAR)
		{
			strain = 0.5f * (F + !F) - id; // Cauchy strain
		}
		else
		{
			strain = 0.5f * (!F * F - id);	// Green strain	
		}

		// compute PK1 stress tensors [Sifakis]
		Matrix3R piolae;
		if (material == MMT_STVK)
		{
			// Green strain - StVK model
			piolae = F * (2 * mu * strain + lambda * strain.Trace() * id);
		}
		else if (material == MMT_LINEAR)
		{
			// Cauchy strain - linear elasticity
			piolae = 2 * mu * strain + lambda * strain.Trace() * id;
		}
		else if (material == MMT_COROTATIONAL)
		{
			piolae = R * (2 * mu * strain + lambda * strain.Trace() * id);
		}
		else if (material == MMT_NEO_HOOKEAN)
		{
			// Neo-Hookean [Sifakis]
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			//if (J < 0)
			//	Printf("negative jacobian :)\n");
			piolae = mu * (F - Finvtr) + lambda * log(J) * Finvtr;
		}
		else if (material == MMT_DISTORTIONAL_NH7)
		{
			// Distortional Neo-Hookean
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			real lambda1 = -2.f / 3.f * mu;
			piolae = mu * (F - Finvtr) + lambda1 * log(J) * Finvtr;
		}
		else if (material == MMT_NEO_HOOKEAN_OGDEN)
		{
			// Neo-Hookean [Smith]->[Ogden]
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			piolae = mu * (F - Finvtr) + lambda * (J - 1) * J * Finvtr;
		}
		else if (material == MMT_DISTORTIONAL_LINEAR)
		{
			// ver 1: zero Poisson ratio
			//piolae = mYoungsModulus * strain;

			// ver 2: linear deviatoric elasticity
			piolae = (2.f * mu) * strain - (2.f * mu / 3.f) * strain.Trace() * id; // deviatoric stress
		}
		else if (material == MMT_DISTORTIONAL_MOONEY)
		{
			// ver 1: Mooney (zero lambda and no J term)
			piolae = mu * F;
		}
		else if (material == MMT_DISTORTIONAL_NH2)
		{
			// ver 2: Mooney with zero Poisson ratio
			piolae = 0.5f * femPhysics->GetYoungsModulus() * F;
		}
		else if (material == MMT_DISTORTIONAL_NH3)
		{
			// ver 3: J-scaled - grows without bound under compression [Smith]->[Bower]
			real J = F.Determinant();
			piolae = mu * pow(J, -2.f / 3.f) * F;
		}
		else if (material == MMT_DISTORTIONAL_BW)
		{
			// ver 4: Bonet & Wood
			real J = F.Determinant();
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			Matrix3R C = !F * F;
			real I1 = C.Trace();
			piolae = mu * pow(J, -2.f / 3.f) * (F - (I1 / 3.f) * Finvtr);
		}
		else if (material == MMT_DISTORTIONAL_OGDEN)
		{
			// ver 5: only the mu terms [Smith]->[Ogden]
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			piolae = mu * (F - Finvtr);
		}
		else if (material == MMT_DISTORTIONAL_NH6)
		{
			// ver 5: only the mu terms of the first material proposed by [Smith]
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			piolae = mu * (F - J * Finvtr);
		}

		return piolae;
	}

	void ElasticEnergy::ComputeForces(const FemPhysicsBase* femPhysics, std::vector<Vector3R>& fout, int material)
	{
		// go through all elements and add up nodal forces
		// TODO: OpenMP acceleration
		for (int e = 0; e < (int)femPhysics->GetNumElements(); e++)
		{
			auto piolae = ComputeElementStress(femPhysics, e, material);
			Matrix3R forces = -femPhysics->GetElementInitialVolume(e) * piolae * femPhysics->GetBarycentricJacobianMatrix(e); // material description as the volume is the undeformed one
			uint32 i[4];
			i[0] = femPhysics->GetGlobalIndex(e, 0);
			i[1] = femPhysics->GetGlobalIndex(e, 1);
			i[2] = femPhysics->GetGlobalIndex(e, 2);
			i[3] = femPhysics->GetGlobalIndex(e, 3);
			for (int j = 1; j < 4; j++)
			{
				Vector3R f = forces(j - 1); // the j-1 column of 'forces'
				fout[i[j]] += f;
				fout[i[0]] -= f;
			}
		}

	}

	void ElasticEnergy::ComputeLocalForceDifferential(const FemPhysicsBase* femPhysics, uint32 e, const Vector3R dx[4], Vector3R df[4])
	{
		// compute deformation gradient increment dF
		const Vector3R& x0 = dx[0];
		const Vector3R& x1 = dx[1];
		const Vector3R& x2 = dx[2];
		const Vector3R& x3 = dx[3];
		Vector3R d1 = x1 - x0;
		Vector3R d2 = x2 - x0;
		Vector3R d3 = x3 - x0;
		Matrix3R mat(d1, d2, d3);
		Matrix3R X = femPhysics->GetBarycentricJacobianMatrix(e);
		Matrix3R dF = mat * !X;

		Matrix3R F;
		femPhysics->ComputeDeformationGradient(e, F); // we don't need it for linear, but never mind

		Matrix3R dP, id;

		// Lame coefficients
		const real mu = femPhysics->GetShearModulus();
		real lambda = femPhysics->GetLameFirstParam();

		if (femPhysics->GetMaterial() == MMT_LINEAR)
		{
			// stress differential - linear elasticity
			dP = mu * (dF + !dF) + lambda * dF.Trace() * id;
		}
		else if (femPhysics->GetMaterial() == MMT_COROTATIONAL)
		{
			// corotational
			Matrix3R R, S;
			ComputePolarDecomposition(F, R, S);
			Matrix3R dS = !R * dF;
			dP = 2 * mu * dF + lambda * dS.Trace() * R;
			// nonlinear correction term - comment below to obtain the implicit corotational method in [Mueller]
			real trS = S.Trace();
			Matrix3R A = S - trS * id;
			Vector3R w(dS(1, 2) - dS(2, 1), dS(0, 2) - dS(2, 0), dS(0, 1) - dS(1, 0));
			Vector3R r = A.GetInverse() * w;
			Matrix3R X = Matrix3R::Skew(r);
			Matrix3R dR = R * X;
			dP = dP + (lambda * (trS - 3) - 2 * mu) * dR;
		}
		else if (femPhysics->GetMaterial() == MMT_STVK)
		{
			Matrix3R E = 0.5f * (!F * F - id);
			Matrix3R dE = 0.5f * (!dF * F + !F * dF);
			dP = dF * (2 * mu * E + lambda * E.Trace() * id)
				+ F * (2 * mu * dE + lambda * dE.Trace() * id);
		}
		else if (femPhysics->GetMaterial() == MMT_NEO_HOOKEAN)
		{
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			dP = mu * dF + (mu - lambda * log(J)) * Finvtr * !dF * Finvtr + lambda * (Finv * dF).Trace() * Finvtr;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_NH7)
		{
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			real lambda1 = -2.f / 3.f * mu;
			dP = mu * dF + (mu - lambda1 * log(J)) * Finvtr * !dF * Finvtr + lambda1 * (Finv * dF).Trace() * Finvtr;
		}
		else if (femPhysics->GetMaterial() == MMT_NEO_HOOKEAN_OGDEN)
		{
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			dP = mu * dF + (mu - lambda * (J - 1)) * Finvtr * !dF * Finvtr + lambda * J * (Finv * dF).Trace() * Finvtr;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_LINEAR)
		{
			// ver 1 : zero Poisson ratio
			//dP = 0.5f * mYoungsModulus * (dF + !dF);

			// ver 2: linear deviatoric elasticity
			dP = mu * (dF + !dF) - (2.f / 3.f) * mu * dF.Trace() * id;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_MOONEY)
		{
			// ver 1: Mooney (zero lambda)
			dP = mu * dF;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_NH2)
		{
			// ver 2: zero Poisson ratio
			dP = 0.5f * femPhysics->GetYoungsModulus() * dF;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_NH3)
		{
			// ver 3: J-scaled
			real J = F.Determinant();
			dP = mu * pow(J, -2.f / 3.f) * dF;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_BW)
		{
			// ver 4: Bonet & Wood
			real J = F.Determinant();			
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			Matrix3R C = !F * F;
			real I1 = C.Trace();
			
			Matrix3R Q = F - (1.f / 3.f) * I1 * Finvtr;
			Matrix3R dC = !dF * F + !F * dF;			
			Matrix3R dQ = dF - (1.f / 3.f) * (dC.Trace() * Finvtr - I1 * Finvtr * !dF * Finvtr);
			dP = mu * pow(J, -2.f / 3.f) * (dQ - (2.f / 3.f) * (Finv * dF).Trace() * Q);
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_OGDEN)
		{
			// only the mu terms
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			dP = mu * dF + mu * Finvtr * !dF * Finvtr;
		}
		else if (femPhysics->GetMaterial() == MMT_DISTORTIONAL_NH6)
		{
			// only the mu terms of the first material proposed by [Smith]
			Matrix3R Finv = F.GetInverse();
			Matrix3R Finvtr = !Finv;
			real J = F.Determinant();
			dP = mu * dF + mu * J * Finvtr * !dF * Finvtr - mu * J * (Finv * dF).Trace() * Finvtr;
		}

		Matrix3R dH = -femPhysics->GetElementInitialVolume(e) * dP * X;

		Vector3R f3;
		for (int j = 1; j < 4; j++)
		{
			Vector3R f = dH(j - 1);
			df[j] = f;
			f3 -= f;
		}
		df[0] = f3;
	}

	void ElasticEnergy::ComputeLocalStiffnessMatrixFromDifferential(const FemPhysicsBase* femPhysics, uint32 e, EigenMatrix& Klocal)
	{
		uint32 numNodes = 4; // local nodes
		uint32 numDofs = numNodes * 3;

		std::vector<Vector3R> df(numNodes);
		std::vector<Vector3R> basisI(numNodes);
		std::vector<real> basisJ(numDofs);
		Klocal.resize(numDofs, numDofs);

		// pre-compute inverse diagonal values
		for (size_t j = 0; j < numDofs; j++)
		{
			// prepare unit vector
			std::fill(basisJ.begin(), basisJ.end(), 0);
			basisJ[j] = 1;
			// compute A-dot product
			Vector3R* dx2 = (Vector3R*)basisJ.data();
			ElasticEnergy::ComputeLocalForceDifferential(femPhysics, e, dx2, &df[0]);
			for (size_t i = 0; i < numDofs; i++)
			{
				// prepare unit vector
				std::fill(basisI.begin(), basisI.end(), Vector3R(0));
				basisI[i / 3][i % 3] = 1;
				real val = InnerProduct<Vector3R, real>(basisI, df);
				Klocal(i, j) = -val;
			}
		}
	}

	void ElasticEnergy::AssembleStiffnessMatrix(const FemPhysicsBase* femPhysics, EigenMatrix& stiffnessMatrix, EigenMatrix* bcStiffnessMatrix)
	{
		PROFILE_SCOPE("Assemble K");

		size_t numNodes = femPhysics->GetNumFreeNodes();
		size_t numDofs = numNodes * 3;
		uint32 numBCs = (uint32)(femPhysics->GetNumNodes() - numNodes);
		stiffnessMatrix.resize(numDofs, numDofs);
		stiffnessMatrix.setZero();
		if (bcStiffnessMatrix)
		{
			bcStiffnessMatrix->resize(numDofs, numBCs * 3);
			bcStiffnessMatrix->setZero();
		}
		// go through all linear elements (tetrahedra)
		for (uint32 i = 0; i < femPhysics->GetNumElements(); i++)
		{
			EigenMatrix Klocal;
			ElasticEnergy::ComputeLocalStiffnessMatrixFromDifferential(femPhysics, i, Klocal);
			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (uint32 j = 0; j < femPhysics->GetNumLocalNodes(); j++)
			{
				uint32 jGlobal = femPhysics->GetGlobalIndex(i, j);
				for (uint32 k = 0; k < femPhysics->GetNumLocalNodes(); k++)
				{
					uint32 kGlobal = femPhysics->GetGlobalIndex(i, k);
					if (jGlobal < numBCs)
						continue;
					int jOffset = (jGlobal - numBCs) * NUM_POS_COMPONENTS;
					int kOffset = (kGlobal - numBCs) * NUM_POS_COMPONENTS;

					// add the the whole 3x3 block to the block matrix
					for (size_t x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (size_t y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							real val = Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
							if (kGlobal < numBCs)
							{
								if (bcStiffnessMatrix)
									bcStiffnessMatrix->coeffRef(jOffset + x, kGlobal * NUM_POS_COMPONENTS + y) += val;
							}
							else
							{
								stiffnessMatrix.coeffRef(jOffset + x, kOffset + y) += val;
							}
						}
					}
				}
			}
		}
	}

	// duplicated code
	void ElasticEnergy::AssembleStiffnessMatrix(const FemPhysicsBase* femPhysics, SparseMatrix& stiffnessMatrix, EigenMatrix* bcStiffnessMatrix)
	{
		PROFILE_SCOPE("Assemble K");

		size_t numNodes = femPhysics->GetNumFreeNodes();
		size_t numDofs = numNodes * 3;
		uint32 numBCs = (uint32)(femPhysics->GetNumNodes() - numNodes);
		stiffnessMatrix.resize(numDofs, numDofs);
		stiffnessMatrix.setZero();
		std::vector<Eigen::Triplet<real>> triplets;
		if (bcStiffnessMatrix)
		{
			bcStiffnessMatrix->resize(numDofs, numBCs * 3);
			bcStiffnessMatrix->setZero();
		}
		// go through all linear elements (tetrahedra)
		for (uint32 i = 0; i < femPhysics->GetNumElements(); i++)
		{
			EigenMatrix Klocal;
			ElasticEnergy::ComputeLocalStiffnessMatrixFromDifferential(femPhysics, i, Klocal);
			// the local stiffness matrix Klocal is organized in 3x3 blocks for each pair of nodes
			for (uint32 j = 0; j < femPhysics->GetNumLocalNodes(); j++)
			{
				uint32 jGlobal = femPhysics->GetGlobalIndex(i, j);
				for (uint32 k = 0; k < femPhysics->GetNumLocalNodes(); k++)
				{
					uint32 kGlobal = femPhysics->GetGlobalIndex(i, k);
					if (jGlobal < numBCs)
						continue;
					int jOffset = (jGlobal - numBCs) * NUM_POS_COMPONENTS;
					int kOffset = (kGlobal - numBCs) * NUM_POS_COMPONENTS;

					// add the the whole 3x3 block to the block matrix
					for (uint32 x = 0; x < NUM_POS_COMPONENTS; x++)
					{
						for (uint32 y = 0; y < NUM_POS_COMPONENTS; y++)
						{
							real val = Klocal(j * NUM_POS_COMPONENTS + x, k * NUM_POS_COMPONENTS + y);
							if (kGlobal < numBCs)
							{
								if (bcStiffnessMatrix)
									bcStiffnessMatrix->coeffRef(jOffset + x, kGlobal * NUM_POS_COMPONENTS + y) += val;
							}
							else
							{
								triplets.push_back(Eigen::Triplet<real>(jOffset + x, kOffset + y, val));
							}
						}
					}
				}
			}
		}
		stiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());
	}

} // namespace FEM_SYSTEM