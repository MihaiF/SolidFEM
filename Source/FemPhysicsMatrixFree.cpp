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

#include "FemPhysicsMatrixFree.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Engine/Profiler.h>
#include "PolarDecomposition.h"
#include "ElasticEnergy.h"
#include "LinearSolver.h"
#include "NewtonSolver.h"
#include "NonlinearCG.h"
#include <iostream>

#pragma warning( disable : 4267) // for size_t-uint conversions

// Bibliography:
// [Sifakis] Sifakis, E., The classical FEM method and discretization methodology, Siggraph course, 2012
// [Teran03] Teran, J. et al., Finite Volume Method for the Simulation of Skeletal Muscle, Eurographics, 2003
// [Teran03] Teran, J. et al., Robust Quasistatic Finite Elements and Flesh Simulation, Eurographics, 2005
// [Mueller] Mueller, M. et al., Real Time Physics Class Notes, Chapter 4 - The finite element method, Siggraph course, 2008
// [Bonet] Bonet, J., Wood, R.D., Nonlinear Continuum Mechanics for Finite Element Analysis
// [Erleben] Erleben, K. et al., Physics-based animation

// We are using explicit FEM: lumped masses (particle masses) and lumped body forces (particle forces)
// Also, it is geometrically linear FEM: all elements are linear displacement/constant strain tetrahedra
// P or PK1 is the first Piola-Kirchoff tensor
// S or PK2 is the second Piola-Kirchoff tensor

//#define CACHED_STIFFNESS_MATRIX

namespace FEM_SYSTEM
{
	// node indices for each face
	int faces[4][3] = {
			{ 1, 2, 3 },
			{ 0, 2, 3 },
			{ 0, 1, 3 },
			{ 0, 1, 2 } };

	FemPhysicsMatrixFree::FemPhysicsMatrixFree(std::vector<Tet>& tetrahedra,
		std::vector<Node>& allNodes, const FemConfig& config)
		: FemPhysicsBase(config)
		, nodes(allNodes)
		, hasCollisions(false)
	{
		tets.resize(tetrahedra.size());
		for (size_t i = 0; i < tets.size(); i++)
		{
			Tetrahedron& tet = tets.at(i);

			tet.i[0] = tetrahedra[i].idx[0];
			tet.i[1] = tetrahedra[i].idx[1];
			tet.i[2] = tetrahedra[i].idx[2];
			tet.i[3] = tetrahedra[i].idx[3];
		}

#ifndef USE_CONSTRAINT_BCS
		ReshuffleFixedNodes();
#else
		mNumBCs = 0;
		mReshuffleMap.resize(nodes.size());
		for (uint32 i = 0; i < nodes.size(); i++)
		{
			mReshuffleMap[i] = i;
			if (nodes[i].invMass == 0)
			{
				nodes[i].invMass = 1;
				AddDirichletBC(i, AXIS_X | AXIS_Y | AXIS_Z);
			}
		}
#endif

		if (config.mCustomConfig != nullptr)
		{
			const Config* cfg = (Config*)config.mCustomConfig;
			mConfig = *cfg;
		}

		// damping params (partial Rayleigh damping: we reuse the matrix structure from the stiffness matrix)
		// TODO: use Rayleigh damping params instead
		ed = config.mDampingYoungsModulus;
		nud = config.mDampingPoissonRatio;
		real omnd = 1.f - nud;
		real om2nd = 1.f - 2 * nud;
		real fd = ed / (1.f + nud) / om2nd;
		Ed = fd * Matrix3R(omnd, nud, nud,
			nud, omnd, nud,
			nud, nud, omnd);

		// prepare the vector of lumped masses
		std::vector<real> masses(nodes.size(), 0.f);

		// init the test and shape matrices
		mTotalInitialVol = 0;
		for (size_t i = 0; i < tets.size(); i++)
		{
			Tetrahedron& tet = tets.at(i);
			const Vector3R& x0 = nodes.at(tet.i[0]).pos0;
			const Vector3R& x1 = nodes.at(tet.i[1]).pos0;
			const Vector3R& x2 = nodes.at(tet.i[2]).pos0;
			const Vector3R& x3 = nodes.at(tet.i[3]).pos0;
			Vector3R d1 = x1 - x0;
			Vector3R d2 = x2 - x0;
			Vector3R d3 = x3 - x0;
			Matrix3R mat(d1, d2, d3); // this is the reference shape matrix Dm [Sifakis][Teran03]
			tet.X = mat.GetInverse(); // Dm^-1
			tet.Xtr = !tet.X; // Dm^-T; pre-stored but could not be for saving memory
			tet.vol = (mat.Determinant()) / 6.f; // signed volume of the tet
			mTotalInitialVol += tet.vol;

			real lumpedMass = 0.25f * tet.vol * config.mDensity;
			masses[tet.i[0]] += lumpedMass;
			masses[tet.i[1]] += lumpedMass;
			masses[tet.i[2]] += lumpedMass;
			masses[tet.i[3]] += lumpedMass;

			// compute face areas
			Vector3R b[4];
			for (int j = 0; j < 4; j++)
			{
				int i1 = faces[j][0];
				int i2 = faces[j][1];
				int i3 = faces[j][2];

				int j1 = tet.i[i1];
				int j2 = tet.i[i2];
				int j3 = tet.i[i3];
				int j4 = tet.i[j];

				tet.NA[j] = cross(nodes[j2].pos - nodes[j1].pos, nodes[j3].pos - nodes[j1].pos);
				if (dot(tet.NA[j], nodes[j4].pos - nodes[j1].pos) > 0)
					tet.NA[j].Flip();

				Vector3R splitFaceForce = (-1.f / 6.f) * tet.NA[j];
				b[i1] += splitFaceForce;
				b[i2] += splitFaceForce;
				b[i3] += splitFaceForce;
			}
			//tet.Bm = Matrix3R(b[1], b[2], b[3]);
			
			// compute the gradient of the shape functions [Erleben]
			Vector3R y[4];
			y[1] = tet.X[0];
			y[2] = tet.X[1];
			y[3] = tet.X[2];
			y[0] = y[1] + y[2] + y[3];
			y[0].Flip();

			// compute the Cauchy strain Jacobian matrix according to [Mueller]: H = de/dx
			for (int j = 0; j < 4; j++)
			{
				tet.Hn[j] = Matrix3R(y[j]);
				tet.Hs[j] = Matrix3R(0, y[j].Z(), y[j].Y(),
					y[j].Z(), 0, y[j].X(),
					y[j].Y(), y[j].X(), 0);
			}
			
			// for linear FEM we can actually precompute the tangent stiffness matrix
			real s = mYoungsModulus / (1.f + mPoissonRatio);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					tet.K[j][k] = tet.vol * (tet.Hn[j] * mNormalElasticityMatrix * tet.Hn[k] + s * tet.Hs[j] * (!tet.Hs[k]));
				}
			}
		}

		// set masses
		for (size_t i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].invMass != 0 && masses[i] != 0)
				nodes[i].invMass = 1.f / masses[i];
		}

		lambdaAcc.resize(tets.size());
		dH.resize(tets.size());

		mForces.resize(nodes.size());

		// build lumped mass matrix
		uint32 numDofs = GetNumFreeNodes() * 3;
		mMassMatrix.resize(numDofs, numDofs);
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			real mass = 1.f / nodes[i + mNumBCs].invMass;
			mMassMatrix.coeffRef(i * 3, i * 3) = mass;
			mMassMatrix.coeffRef(i * 3 + 1, i * 3 + 1) = mass;
			mMassMatrix.coeffRef(i * 3 + 2, i * 3 + 2) = mass;
		}
	}

	void FemPhysicsMatrixFree::ReshuffleFixedNodes()
	{
		std::vector<Node> newNodes;
		mReshuffleMap.resize(nodes.size()); // map from old indices to new ones
		// add fixed nodes first
		mNumBCs = 0;
		for (size_t i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].invMass == 0)
			{
				mReshuffleMap[i] = newNodes.size();
				newNodes.push_back(nodes[i]);
				mNumBCs++;
			}
		}
		// then the other nodes
		for (size_t i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].invMass != 0)
			{
				mReshuffleMap[i] = newNodes.size();
				newNodes.push_back(nodes[i]);
			}
		}
		nodes = newNodes; // replace old nodes with shuffled ones
		// remap tets
		originalTets = tets;
		for (size_t i = 0; i < tets.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tets[i].i[j] = mReshuffleMap[tets[i].i[j]];
			}
		}
	}

	void FemPhysicsMatrixFree::Step(real dt)
	{
		if (mSimType == ST_IMPLICIT)
		{
			// save positions before step - breaks dynamic Dirichlet BCs!
			for (size_t i = 0; i < GetNumNodes(); i++)
			{
				nodes[i].pos0 = nodes[i].pos;
			}
			
			real h = dt / mNumSteps;
			for (int i = 0; i < mNumSteps; i++)
			{
				mTimeStep = h;
				mForceFraction = 1;
				Solve();
				HandleCollisions(h);
			}
		}
		else if (mSimType == ST_EXPLICIT)
		{
			real h = dt / mNumSteps;
			for (int i = 0; i < mNumSteps; i++)
			{
				SubStep(h);
			}
		}
		else if (mSimType == ST_STATIC)
		{
			if (mForceFraction == 0)
			{
				mForceFraction = 1;
				Solve();				
			}
		}
		else if (mSimType == ST_QUASI_STATIC)
		{
			if (mForceFraction < 1)
			{	
				mForceFraction = std::min(real(1), mForceFraction + mForceStep);
				Solve();
			}
		}

	}

	void FemPhysicsMatrixFree::UpdatePositions(std::vector<Node>& newNodes)
	{
		// TODO: use mNumBCs and mReshuffleMapInv
		for (size_t i = 0; i < newNodes.size(); i++)
		{
			if (newNodes[i].invMass == 0)
				nodes[mReshuffleMap[i]].pos = newNodes[i].pos;
		}
	}

	// Helper function to compute the nodal elastic forces given the current configuration of the system
	// It relies on linear elasticity (small strain tensor) to precompute the tangent stiffness matrix [Mueller]
	void FemPhysicsMatrixFree::ComputeForces_LinearElasticity()
	{
		// TODO: Rayleigh damping D = alpha * M + beta * K (does not work for explicit integration at least)
		// TODO: Co-rotational formulation using the pre-stored element stiffness matrix

		disp.resize(nodes.size());
		for (size_t i = 0; i < nodes.size(); i++)
		{
			disp[i] = nodes[i].pos - nodes[i].pos0; // + beta * nodes[i].vel;
		}

		for (size_t i = 0; i < tets.size(); i++)
		{
			Tetrahedron& tet = tets.at(i);

			// compute node forces
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					nodes[tet.i[j]].force -= tet.K[j][k] * disp[tet.i[k]];
				}
				//if (nodes[nod].invMass != 0)
					//nodes[nod].force -= (alpha / nodes[nod].invMass) * nodes[nod].vel; // Rayleigh damping: D = alpha * M (why the minus?)
			}
		}
	}

	// Helper function to compute the nodal elastic forces for a given tetrahedron and rotation
	// It uses a co-rotated linear elasticity (a.k.a. warped stiffness) formulation that precomputes the tangent stiffness matrix [Mueller]
	void FemPhysicsMatrixFree::ComputeForces_CorotationalElasticity(const Tetrahedron& tet, const Matrix3R& R)
	{
		// compute internal forces
		Vector3R u[4];
		const real alpha = 0.2f;
		const real beta = 0.016f;
		// compute displacements
		for (int j = 0; j < 4; j++)
		{
			int n = tet.i[j];
			u[j] = !R * nodes.at(n).pos - nodes.at(n).pos0 + beta * nodes.at(n).vel; // Rayleigh damping: D = beta * K; TODO: precompute
		}

		// compute node forces
		for (int j = 0; j < 4; j++)
		{
			Vector3R f;
			for (int k = 0; k < 4; k++)
			{
				f -= tet.K[j][k] * u[k];
			}
			int nod = tet.i[j];
			if (nodes.at(nod).invMass != 0)
				f -= (alpha / nodes.at(nod).invMass) * nodes.at(nod).vel; // Rayleigh damping: D = alpha * M (why the minus?)
			nodes.at(nod).force += R * f;
		}
	}

	void FemPhysicsMatrixFree::ComputeDeformationGradient(uint32 e, Matrix3R& F) const
	{
		const Tetrahedron& tet = tets[e];
		// compute deformed/spatial shape matrix Ds [Sifakis]
		const Vector3R& x0 = nodes.at(tet.i[0]).pos;
		const Vector3R& x1 = nodes.at(tet.i[1]).pos;
		const Vector3R& x2 = nodes.at(tet.i[2]).pos;
		const Vector3R& x3 = nodes.at(tet.i[3]).pos;
		Vector3R d1 = x1 - x0;
		Vector3R d2 = x2 - x0;
		Vector3R d3 = x3 - x0;
		Matrix3R Ds(d1, d2, d3);
		// compute deformation gradient
		F = Ds * tet.X;
	}

	// Helper function to compute the nodal elastic forces given the current configuration of the system
	// STVK = St. Vennant - Kirchoff material model - just use Green strain instead of small strain, but keep linear elasticity law
	// TODO: remove as obsolete
	template<int NodalForces>
	void FemPhysicsMatrixFree::ComputeForces()
	{
		// compute internal forces: reset nodal forces and add gravity to non-fixed nodes
		for (size_t i = 0; i < nodes.size(); i++)
		{
			nodes.at(i).force.SetZero();
		}
		
#ifdef CACHED_STIFFNESS_MATRIX
		if (Material == MMT_LINEAR)
		{
			// a faster way of computing linear elastic nodal forces: f = Ku (K is cached)
			ComputeForces_LinearElasticity();
			return; // we bypass the rest of the code
		}
#endif

		// helper constants
		const Matrix3R id;
		const real t = ed / (1.f + nud);
		const real mu = GetShearModulus();
		const real lambda = GetLameFirstParam();
		const real mud = 0.5f * t;
		const real lambdad = t * nud / (1.f - 2 * nud);

		// go through all elements and add up nodal forces
		for (size_t i = 0; i < tets.size(); i++)
		{
			// compute nodal forces for this element only
			Tetrahedron& tet = tets.at(i);

			Matrix3R F;
			ComputeDeformationGradient(i, F);

			Matrix3R R, U;
			if (mMaterial == MMT_COROTATIONAL)
			{
				ComputePolarDecomposition(F, R, U);
				F = U; // this seems to be broken for the torus!

#ifdef CACHED_STIFFNESS_MATRIX
				// precomputed K matrix alternative
				ComputeForces_CorotationalElasticity(tet, R);
				continue;
#endif
			}

			// compute V = dF/dt - the velocity gradient
			const Vector3R& v0 = nodes.at(tet.i[0]).vel;
			const Vector3R& v1 = nodes.at(tet.i[1]).vel;
			const Vector3R& v2 = nodes.at(tet.i[2]).vel;
			const Vector3R& v3 = nodes.at(tet.i[3]).vel;
			Matrix3R V = Matrix3R(v1 - v0, v2 - v0, v3 - v0) * tet.X;

			// compute strain and strain rate as tensors
			Matrix3R strain, strainRate;
			if (mMaterial == MMT_LINEAR || mMaterial == MMT_COROTATIONAL || mMaterial == MMT_DISTORTIONAL_LINEAR)
			{
				strainRate = 0.5f * (V + !V); // Cauchy strain rate
				strain = 0.5f * (F + !F) - id; // Cauchy strain
			}
			else
			{
				strain = 0.5f * (!F * F - id);	// Green strain	
				strainRate = 0.5f * (!F * V + !V * F); // Green strain rate
			}

			// [Sifakis] version of computing the stress tensor
			// compute PK1 stress tensors
			Matrix3R piolae, piolad;
			if (mMaterial == MMT_STVK)
			{
				// Green strain - StVK model
				piolae = F * (2 * mu * strain + lambda * strain.Trace() * id);
				piolad = F * (mud * strainRate + 0.5f * lambdad * strainRate.Trace() * id);
			}
			else if (mMaterial == MMT_LINEAR || mMaterial == MMT_DISTORTIONAL_LINEAR)
			{
				// Cauchy strain - linear elasticity
				piolae = 2 * mu * strain + lambda * strain.Trace() * id;
				piolad = mud * strainRate + 0.5f * lambdad * strainRate.Trace() * id;
			}
			else if (mMaterial == MMT_COROTATIONAL)
			{
				piolae = R * (2 * mu * strain + lambda * strain.Trace() * id);
				piolad = R * (mud * strainRate + 0.5f * lambdad * strainRate.Trace() * id);
			}
			else if (mMaterial == MMT_NEO_HOOKEAN)
			{
				// Neo-Hookean [Sifakis]
				Matrix3R Finv = F.GetInverse();
				Matrix3R Finvtr = !Finv;
				real J = F.Determinant();
				piolae = mu * (F - Finvtr) + lambda * log(J) * Finvtr;
				piolad = Matrix3R::Zero();
			}

			// compute the total PK1 stress tensor
			Matrix3R P = piolae + piolad;

			if (NodalForces == NFT_DIRECT)
			{
				// [Sifakis] approach at computing nodal forces 
				// P needs to be the first Piola-Kirchoff (PK1) stress tensor
				Matrix3R forces = -tet.vol * P * tet.Xtr; // material description as the volume is the undeformed one
				//Matrix3R forces = P * tet.Bm; // the matrix from [Teran03], it should be Bm = -V * Dm^{-T} [Teran05]
				for (int j = 1; j < 4; j++)
				{
					Vector3R f = forces(j - 1); // the j-1 column of 'forces'
					nodes.at(tet.i[j]).force += f;
					nodes.at(tet.i[0]).force -= f;
				}
			}
			else
			{
				// finite volume approach to nodal forces [Mueller][Teran03]
				for (int j = 0; j < 4; j++)
				{
					int j1 = tet.i[faces[j][0]];
					int j2 = tet.i[faces[j][1]];
					int j3 = tet.i[faces[j][2]];
					int j4 = tet.i[j];

					// when using deformed area the stress tensor should be Cauchy stress (spatial setting) [Mueller]
					//Vector3R normalArea = cross(nodes[j2].pos - nodes[j1].pos, nodes[j3].pos - nodes[j1].pos);
					//if (dot(normalArea, nodes[j4].pos - nodes[j1].pos) > 0) // if pointing in the same direction
					//	normalArea.Flip();
					//Vector3R splitFaceForce = (-1.f / 6.f) * P * normalArea;

					// there is a typo in [Mueller]: it should be 6 instead of 3 in the denominator
					// see [Teran03] just before sub-section 4.1
					Vector3R splitFaceForce = (-1.f / 6.f) * P * tet.NA[j];
					nodes[j1].force += splitFaceForce;
					nodes[j2].force += splitFaceForce;
					nodes[j3].force += splitFaceForce;
				}

				// alternatively we can precompute the nodal vectors, but that amounts to the Bm matrix approach above
				//for (int j = 1; j < 4; j++)
				//{
				//	Vector3R force = P * tet.b[j];
				//	nodes[tet.i[j]].force += force;
				//	nodes[tet.i[0]].force -= force;
				//}
			}
		}
	}

	// explicit instantiations for outside use
	template void FemPhysicsMatrixFree::ComputeForces<FemPhysicsMatrixFree::NFT_DIRECT>();
	template void FemPhysicsMatrixFree::ComputeForces<FemPhysicsMatrixFree::NFT_FINITE_VOLUME>();

	// simplest way of doing time dependent FEM - explicit integration of the discretized system
	void FemPhysicsMatrixFree::SubStep(real h)
	{
		// accumulate gravity and elastic forces in the nodes
		ComputeForces<NFT_DIRECT>();

		// integrate node velocities and positions using Symplectic Euler
		for (size_t i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].invMass == 0)
				continue;
			nodes.at(i).vel += (h * nodes.at(i).invMass) * nodes.at(i).force + h * mGravity;
			nodes.at(i).pos += h * nodes.at(i).vel;
		}
	}

	inline void AddMatrix3(EigenMatrix& K, int n, int m, const Matrix3R& Ke, real h2)
	{
		int nn = 3 * n;
		int mm = 3 * m;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
#ifdef SPARSE
				K.coeffRef(nn + i, mm + j) -= h2 * Ke.m[i][j];
#else
				K(nn + i, mm + j) -= h2 * Ke.m[i][j];
#endif
			}
		}
	}

	void FemPhysicsMatrixFree::ComputeForceDifferential(const std::vector<Vector3R>& dx, std::vector<Vector3R>& df) const
	{
		PROFILE_SCOPE("Differential");
		memset(&df[0], 0, df.size() * sizeof(Vector3R));

		// Lame coefficients
		const real mu = GetShearModulus();
		const real lambda = GetLameFirstParam();

		Matrix3R id;
		Vector3R zero;
		#pragma omp parallel for
		for (int i = 0; i < (int)tets.size(); i++)
		{
			const Tetrahedron& tet = tets.at(i);
			// compute deformation gradient increment dF
			const Vector3R& x0 = tet.i[0] < mNumBCs ? zero : dx[tet.i[0] - mNumBCs];
			const Vector3R& x1 = tet.i[1] < mNumBCs ? zero : dx[tet.i[1] - mNumBCs];
			const Vector3R& x2 = tet.i[2] < mNumBCs ? zero : dx[tet.i[2] - mNumBCs];
			const Vector3R& x3 = tet.i[3] < mNumBCs ? zero : dx[tet.i[3] - mNumBCs];
			Vector3R d1 = x1 - x0;
			Vector3R d2 = x2 - x0;
			Vector3R d3 = x3 - x0;
			Matrix3R mat(d1, d2, d3);
			Matrix3R dF = mat * tet.X;

			Matrix3R F;
			ComputeDeformationGradient(i, F); // we don't need it for linear, but never mind

			Matrix3R dP;
			if (mMaterial == MMT_LINEAR || mMaterial == MMT_DISTORTIONAL_LINEAR)
			{
				// stress differential - linear elasticity
				dP = mu * (dF + !dF) + lambda * dF.Trace() * id;
			}
			else if (mMaterial == MMT_COROTATIONAL)
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
			else if (mMaterial == MMT_STVK)
			{
				Matrix3R E = 0.5f * (!F * F - id);
				Matrix3R dE = 0.5f * (!dF * F + !F * dF);
				dP = dF * (2 * mu * E + lambda * E.Trace() * id)
					+ F * (2 * mu * dE + lambda * dE.Trace() * id);
			}
			else if (mMaterial == MMT_NEO_HOOKEAN)
			{
				Matrix3R Finv = F.GetInverse();
				Matrix3R Finvtr = !Finv;
				real J = F.Determinant();
				dP = mu * dF + (mu - lambda * log(J)) * Finvtr * !dF * Finvtr + lambda * (Finv * dF).Trace() * Finvtr;
			}

			// compute the double contraction of dF and dP to check PD
			//real indicator = Matrix3R::DoubleContraction(dF, dP);
			//if (indicator < 0)
			//{
			//	Printf("Stiffness matrix is not PD (dF:dP = %g).\n", indicator);
			//}

			dH[i] = -tet.vol * dP * tet.Xtr;
		}

		for (size_t i = 0; i < tets.size(); i++)
		{
			const Tetrahedron& tet = tets.at(i);
			Vector3R f3;
			for (int j = 1; j < 4; j++)
			{
				Vector3R f = dH[i](j - 1);
				if (tet.i[j] >= mNumBCs)
					df[tet.i[j] - mNumBCs] += f;
				f3 -= f;
			}
			if (tet.i[0] >= mNumBCs)
				df[tet.i[0]- mNumBCs] += f3;
		}
	}

	real FemPhysicsMatrixFree::GetTotalVolume() const
	{
		real totalVol = 0;
		for (size_t e = 0; e < tets.size(); e++)
		{
			// compute current volume
			const Tetrahedron& tet = tets.at(e);
			int i0 = tet.i[0];
			int i1 = tet.i[1];
			int i2 = tet.i[2];
			int i3 = tet.i[3];
			const Vector3R& x0 = nodes.at(i0).pos;
			const Vector3R& x1 = nodes.at(i1).pos;
			const Vector3R& x2 = nodes.at(i2).pos;
			const Vector3R& x3 = nodes.at(i3).pos;
			Vector3R d1 = x1 - x0;
			Vector3R d2 = x2 - x0;
			Vector3R d3 = x3 - x0;
			Matrix3R mat(d1, d2, d3); // this is the spatial shape matrix Dm [Sifakis][Teran03]
			real vol = (mat.Determinant()) / 6.f; // signed volume of the tet
			totalVol += vol;
		}
		return totalVol;
	}

	void FemPhysicsMatrixFree::MatrixVectorMultiply(const std::vector<Vector3R>& d, std::vector<Vector3R>& df) const
	{
		PROFILE_SCOPE("Matrix mult");
		// compute M * d / h^2 - K * d (only the second term for quasi-static, identified by mTimeStep = 0)
		ComputeForceDifferential(d, df);
		const real invHSqr = mTimeStep == 0 ? 0 : 1.f / (mTimeStep * mTimeStep);
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
		{
			real mass = 1.f / nodes[i + mNumBCs].invMass;
			df[i] = invHSqr * mass * d[i] - df[i];
		}
	}

	void FemPhysicsMatrixFree::ComputeGradients(std::vector<Vector3R>& r)
	{
		PROFILE_SCOPE("Gradients");
		// compute negative gradient r = f(x, v) + M * g
		for (uint32 i = 0; i < mForces.size(); i++)
			mForces[i].SetZero();
		ComputePressureForces(mForces, mTractionStiffnessMatrix);
		ElasticEnergy::ComputeForces(this, mForces);
		const real invH = mTimeStep == 0 ? 0 : 1.f / mTimeStep;
		const real invHSqr = mTimeStep == 0 ? 0 : 1.f / (mTimeStep * mTimeStep);
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
		{
			if (nodes[i + mNumBCs].invMass == 0)
				continue; // this is for USE_CONSTRAINT_BCS but does no harm
			real mass = 1.f / nodes[i + mNumBCs].invMass;
			r[i] = mForces[i + mNumBCs] + mass * (mForceFraction * mGravity + invH * nodes[i + mNumBCs].vel);
			r[i] -= mass * invHSqr * (nodes[i + mNumBCs].pos - nodes[i + mNumBCs].pos0); // the inertial term
		}
	}

	void FemPhysicsMatrixFree::UpdatePosAndComputeGradients(const std::vector<Vector3R>& pos, std::vector<Vector3R>& r)
	{
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
			nodes[i + mNumBCs].pos = pos[i];

		ComputeGradients(r);
	}

	bool FemPhysicsMatrixFree::Solve()
	{
		MEASURE_TIME("Nonlinear FEM solve");
		PROFILE_SCOPE("Solve");

		bool ret = true;
		if (mConfig.mSolver == NST_NEWTON)
			SolveNewton();
		else if (mConfig.mSolver == NST_NEWTON_LS)
			SolveNewtonLS();
		else if (mConfig.mSolver == NST_NEWTON_CG)
			SolveNewtonCG();
		else if (mConfig.mSolver == NST_NONLINEAR_CG)
			SolveNonlinearConjugateGradient();
		else if (mConfig.mSolver == NST_STEEPEST_DESCENT)
			SolveNonlinearSteepestDescent();
		else if (mConfig.mSolver == NST_GRADIENT_DESCENT)
			SolveGradientDescent(mConfig.mDescentRate);

		if (mTimeStep != 0 && mSimType == ST_IMPLICIT)
		{
			// compute new velocities
			const real invH = 1.f / mTimeStep;
			for (size_t i = 0; i < GetNumFreeNodes(); i++)
			{
				nodes.at(i + mNumBCs).vel = invH * (nodes[i + mNumBCs].pos - nodes[i + mNumBCs].pos0);
			}
		}

		//Printf("final energy: %g\n", ComputeEnergy());

		return ret;
	}

	void FemPhysicsMatrixFree::SolveNewtonCG()
	{
		std::vector<Vector3R> delta(GetNumFreeNodes()); // allocation!
		std::vector<Vector3R> b(GetNumFreeNodes()); // allocation!

		// no need for Newton for linear elasticity
		int numIters = /*mMaterial == MMT_LINEAR ? 1 : */mOuterIterations;

		// Newton steps (outer loop)
		for (int iter = 0; iter < numIters; iter++)
		{
			ComputeGradients(b);

			// solve K * dx = b by Conjugate Gradient
			std::fill(delta.begin(), delta.end(), Vector3R()); // reset the guess to zero; TODO: warm start
			SolveConjugateGradientMF<real, Vector3R, FemPhysicsMatrixFree>(*this, b, delta, mInnerIterations);

			// add dx
			for (size_t i = 0; i < GetNumFreeNodes(); i++)
			{
				//mTotalDisplacements[i] += delta[i];
				nodes.at(i + mNumBCs).pos += delta[i];
			}
		}
	}

	real FemPhysicsMatrixFree::MeritResidual(const EigenVector& rhs)
	{
		if (mConfig.mOptimizer)
			return ComputeEnergy();
		else
			return 0.5 * rhs.squaredNorm();
	}

	EigenVector FemPhysicsMatrixFree::ComputeRhs(const EigenVector& sol)
	{
		// set node positions
		auto solVecs = GetStdVector(sol);
		ASSERT(solVecs.size() == GetNumFreeNodes());
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			nodes[i + mNumBCs].pos = solVecs[i];
		}

		std::vector<Vector3R> b(GetNumFreeNodes()); // allocation!
		ComputeGradients(b);
		EigenVector rhs = GetEigenVector(b);

		// add penalty Dirichlet BC forces
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

			EigenVector dirichletErrors(count);
			count = 0;
			for (uint32 i = 0; i < mDirichletIndices.size(); i++)
			{
				// the i'th BC				
				uint32 flags = mDirichletAxes[i];
				uint32 idx0 = mDirichletIndices[i];
				uint32 idx = idx0 - mNumBCs; // affecting node idx
				// compute the BC error
				if (flags & AXIS_X)
				{
					dirichletErrors(count++) = nodes[idx].pos.x - nodes[idx].pos0.x;
				}
				if (flags & AXIS_Y)
				{
					dirichletErrors(count++) = nodes[idx].pos.y - nodes[idx].pos0.y;
				}
				if (flags & AXIS_Z)
				{
					dirichletErrors(count++) = nodes[idx].pos.z - nodes[idx].pos0.z;
				}
			}
			//Printf("Dirichlet BCs error: %g\n", dirichletErrors.lpNorm<Eigen::Infinity>());
			rhs -= mDirichletStiffness * mDirichletJacobian.transpose() * dirichletErrors;
			// TODO: we can build the rhs term directly without the Jacobian
		}

		//Printf("E=%g\n", ComputeEnergy());
		return rhs;
	}

	template<class MATRIX>
	void FemPhysicsMatrixFree::ComputeSystemMatrix(const EigenVector& s, const EigenVector& y, MATRIX& K)
	{		
		ElasticEnergy::AssembleStiffnessMatrix(this, K);
		//K -= mTractionStiffnessMatrix;
		if (mTimeStep != 0 && mSimType == ST_IMPLICIT)
			K += (1.f / mTimeStep / mTimeStep) * mMassMatrix;

		if (!mDirichletIndices.empty())
		{
			K += mDirichletStiffness * mDirichletJacobian.transpose() * mDirichletJacobian;
		}
		if (mConfig.mOptimizer)
			Printf("energy: %g\n", ComputeEnergy());
	}

	EigenVector FemPhysicsMatrixFree::SolveLinearSystem(EigenMatrix& K, const EigenVector& rhs, const EigenVector& s, const EigenVector& y)
	{
		LinearSolver solver;
		solver.Init(K, LST_LU_PARDISO);
		//solver.SetTolerance(0.1);
		//solver.Init(K, LST_CG);
		return solver.Solve(rhs);
	}

	EigenVector FemPhysicsMatrixFree::SolveLinearSystem(SparseMatrix& K, const EigenVector& rhs, const EigenVector& s, const EigenVector& y)
	{
		LinearSolver solver;
		solver.Init(K, LST_LU_PARDISO);
		return solver.Solve(rhs);
	}

	bool FemPhysicsMatrixFree::SolveNewtonLS()
	{
		AssembleDynamicContributions();
		
		NewtonSolverBackTrack<FemPhysicsMatrixFree, EigenMatrix> solver;
		solver.mNumIterations = mOuterIterations;
		solver.mVerbose = mVerbose ? VL_MINIMUM : VL_NONE;
		solver.mResidualThreshold = mAbsNewtonResidualThreshold;
		solver.mUseProblemSolver = true;
		if (mConfig.mOptimizer)
		{
			solver.mLSCondition = LSC_ARMIJO;
			//solver.mAlpha = 1e-3;
		}
		//solver.mUseBFGS = false;

		// prepare the initial guess with the current configuration
		uint32 size = GetNumFreeNodes() * NUM_POS_COMPONENTS;
		EigenVector solution(size);
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			const Vector3R& p = nodes[i + mNumBCs].pos;
			solution(i * NUM_POS_COMPONENTS) = p.x;
			solution(i * NUM_POS_COMPONENTS + 1) = p.y;
			solution(i * NUM_POS_COMPONENTS + 2) = p.z;
		}

		solver.Solve(*this, size, solution);
		if (mVerbose)
		{
			real vol = GetTotalVolume();
			Printf("vol err: %.2f%%\n", abs(vol - mTotalInitialVol) / mTotalInitialVol * 100);
		}

		// set node positions
		auto solVecs = GetStdVector(solution);
		ASSERT(solVecs.size() == GetNumFreeNodes());
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			nodes[i + mNumBCs].pos = solVecs[i];
		}

		CheckForInversion(true);

		return true;
	}

	bool FemPhysicsMatrixFree::SolveNewton()
	{
		AssembleDynamicContributions();

		NewtonSolver<FemPhysicsMatrixFree, SparseMatrix> solver;
		solver.mNumIterations = mOuterIterations;
		solver.mVerbose = mVerbose ? VL_MINIMUM : VL_NONE;
		solver.mResidualThreshold = mAbsNewtonResidualThreshold;
		solver.mSolverType = LST_LU_PARDISO;
		//solver.mUseFiniteDiff = true;

		// prepare the initial guess with the current configuration
		uint32 size = GetNumFreeNodes() * NUM_POS_COMPONENTS;
		//Printf("#dofs: %d\n", size);
		EigenVector solution(size);
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			const Vector3R& p = nodes[i + mNumBCs].pos;
			solution(i * NUM_POS_COMPONENTS) = p.x;
			solution(i * NUM_POS_COMPONENTS + 1) = p.y;
			solution(i * NUM_POS_COMPONENTS + 2) = p.z;
		}

		solver.Solve(*this, size, solution);
		if (mVerbose)
		{
			real vol = GetTotalVolume();
			Printf("vol err: %.2f%%\n", abs(vol - mTotalInitialVol) / mTotalInitialVol * 100);
		}

		// set node positions
		auto solVecs = GetStdVector(solution);
		ASSERT(solVecs.size() == GetNumFreeNodes());
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
		{
			nodes[i + mNumBCs].pos = solVecs[i];
		}

		CheckForInversion(true);

		return true;
	}

	void FemPhysicsMatrixFree::SolveGradientDescent(real alpha)
	{
		std::vector<Vector3R> b(GetNumFreeNodes()); // allocation!
		for (uint32 iter = 0; iter < mOuterIterations; iter++)
		{
			ComputeGradients(b);

			// perform the gradient descent step (very similar to explicit integration)
			for (size_t i = 0; i < GetNumFreeNodes(); i++)
			{
				Vector3R delta = alpha * b[i];
				nodes[i + mNumBCs].pos += delta;
			}
		}
	}

	void FemPhysicsMatrixFree::SolveNonlinearSteepestDescent()
	{
		// the residual vector
		std::vector<Vector3R> r(GetNumFreeNodes()); // allocation!
		// the differential vector
		std::vector<Vector3R> df(GetNumFreeNodes()); // allocation!

		for (uint32 iter = 0; iter < mOuterIterations; iter++)
		{
			ComputeGradients(r);
			MatrixVectorMultiply(r, df);
			real delta = InnerProduct<Vector3R, real>(r, r);
			if (delta == 0)
				break;
			real alpha = delta / InnerProduct<Vector3R, real>(r, df);
			if (alpha == 0)
				break;

			// perform the gradient descent step (very similar to explicit integration)
			for (size_t i = 0; i < GetNumFreeNodes(); i++)
			{
				Vector3R delta = alpha * r[i];
				nodes[i + mNumBCs].pos += delta;
			}
		}
	}

	void FemPhysicsMatrixFree::SolveNonlinearConjugateGradient()
	{
		std::vector<Vector3R> pos(GetNumFreeNodes()); // allocation!
		for (uint32 i = 0; i < GetNumFreeNodes(); i++)
			pos[i] = nodes[i + mNumBCs].pos;

		NonlinearConjugateGradientMinimizer<FemPhysicsMatrixFree, std::vector<Vector3R>> minimizer;
		minimizer.mOuterIterations = mOuterIterations;
		minimizer.mInnerIterations = mInnerIterations;
		minimizer.mAbsResidualThreshold = mAbsNewtonResidualThreshold;
		minimizer.Solve(*this, GetNumFreeNodes(), pos);
	}

	real FemPhysicsMatrixFree::ComputeEnergy() const
	{
		real val = ElasticEnergy::ComputeEnergy(this);
		// add the gravitational and kinetic part
		real invHSqr = mTimeStep != 0 ? 0.5f / mTimeStep / mTimeStep : 0;
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
		{
			real mass = 1.f / nodes[i + mNumBCs].invMass;
			val -= mass * mForceFraction * mGravity.y * nodes[i + mNumBCs].pos.y; // gravitational
			val -= invHSqr * mass * GetTotalDisplacement(i).LengthSquared(); // kinetic
		}
		// add the pressure part (WIP)
		for (size_t i = 0; i < mTractionSurface.size() / 3; i += 3)
		{
			int base = i * 3;
			// global (shuffled) indices of nodes
			uint32 i1 = mTractionSurface[base];
			uint32 i2 = mTractionSurface[base + 1];
			uint32 i3 = mTractionSurface[base + 2];

			// compute triangle area and normal
			const Vector3R& p1 = GetDeformedPosition(i1);
			const Vector3R& p2 = GetDeformedPosition(i2);
			const Vector3R& p3 = GetDeformedPosition(i3);

			// compute volume formed with the origin
			real vol = (1.f / 6.f) * triple(p1, p2, p3);
			val += abs(vol) * mAppliedPressure;
		}

		//Printf("energy: %g\n", val);
		return val;
	}

	real FemPhysicsMatrixFree::ComputeEnergy(int level) const
	{
		// level 0 -> elastic
		// level 1 -> gravitational
		// level 2 -> vol
		// level 3 -> kinetic

		real val = 0;
		
		if (level == 0)
			val += ElasticEnergy::ComputeEnergy(this);
		
		// add the gravitational and kinetic part
		real invHSqr = mTimeStep != 0 ? 0.5f / mTimeStep / mTimeStep : 0;
		
		for (size_t i = 0; i < GetNumFreeNodes(); i++)
		{
			real mass = 1.f / nodes[i + mNumBCs].invMass;
			if (level == 1)
				val += mass * mForceFraction * mGravity.y * nodes[i + mNumBCs].pos.y; // gravitational
			
			if (level == 3)
				val += invHSqr * mass * GetTotalDisplacement(i).LengthSquared(); // kinetic
		}


		//// add the pressure part (WIP)
		//for (size_t i = 0; i < mTractionSurface.size() / 3; i += 3)
		//{
		//	int base = i * 3;
		//	// global (shuffled) indices of nodes
		//	uint32 i1 = mTractionSurface[base];
		//	uint32 i2 = mTractionSurface[base + 1];
		//	uint32 i3 = mTractionSurface[base + 2];

		//	// compute triangle area and normal
		//	const Vector3R& p1 = GetDeformedPosition(i1);
		//	const Vector3R& p2 = GetDeformedPosition(i2);
		//	const Vector3R& p3 = GetDeformedPosition(i3);

		//	// compute volume formed with the origin
		//	real vol = (1.f / 6.f) * triple(p1, p2, p3);
		//	val += abs(vol) * mAppliedPressure;
		//}

		//Printf("energy: %f\n", val);
		return val;
	}

	void FemPhysicsMatrixFree::ComputePressureForces(Vector3Array& fout, EigenMatrix& Kout) const
	{
		if (mTractionSurface.empty())
			return;

		if (mUseImplicitPressureForces)
		{
			uint32 numDofs = GetNumFreeNodes() * NUM_POS_COMPONENTS;
			Kout.resize(numDofs, numDofs);
			Kout.setZero();
		}

		auto computeForce = [&](const Vector3R& p1, const Vector3R& p2, const Vector3R& p3)->Vector3R
		{
			Vector3R a = p2 - p1;
			Vector3R b = p3 - p1;
			Vector3R normal = cross(a, b);
			real area = 0.5f * normal.Length();
			normal.Normalize();

			Vector3R traction = mAppliedPressure * normal;
			Vector3R force = (area / 3.0f) * traction;
			return force;
		};

		// go through all inner boundary triangles
		for (uint32 t = 0; t < mTractionSurface.size() / 3; t++)
		{
			int base = t * 3;
			// global (shuffled) indices of nodes
			uint32 i1 = mTractionSurface[base];
			uint32 i2 = mTractionSurface[base + 1];
			uint32 i3 = mTractionSurface[base + 2];

			// compute triangle area and normal
			const Vector3R& p1 = GetDeformedPosition(i1);
			const Vector3R& p2 = GetDeformedPosition(i2);
			const Vector3R& p3 = GetDeformedPosition(i3);

			// apply traction to triangle nodes
			Vector3R force = mForceFraction * computeForce(p1, p2, p3);
			i1 = mReshuffleMap[i1];
			i2 = mReshuffleMap[i2];
			i3 = mReshuffleMap[i3];
			fout[i1] += force;
			fout[i2] += force;
			fout[i3] += force;

			if (mUseImplicitPressureForces)
			{
				// compute local stiffness matrix
				Matrix3R K[3];
				Vector3R a = p2 - p1;
				Vector3R b = p3 - p1;
				a.Scale(mAppliedPressure / 6.0f);
				b.Scale(mAppliedPressure / 6.0f);
				K[2] = Matrix3R::Skew(a);
				K[1] = Matrix3R::Skew(-b);
				K[0] = -K[1] - K[2];

				i1 -= mNumBCs;
				i2 -= mNumBCs;
				i3 -= mNumBCs;

				// assemble global stiffness matrix
				for (uint32 x = 0; x < 3; x++)
				{
					for (uint32 y = 0; y < 3; y++)
					{
						Kout(i1 * 3 + x, i1 * 3 + y) += K[0](x, y);
						Kout(i1 * 3 + x, i2 * 3 + y) += K[1](x, y);
						Kout(i1 * 3 + x, i3 * 3 + y) += K[2](x, y);
						Kout(i2 * 3 + x, i1 * 3 + y) += K[0](x, y);
						Kout(i2 * 3 + x, i2 * 3 + y) += K[1](x, y);
						Kout(i2 * 3 + x, i3 * 3 + y) += K[2](x, y);
						Kout(i3 * 3 + x, i1 * 3 + y) += K[0](x, y);
						Kout(i3 * 3 + x, i2 * 3 + y) += K[1](x, y);
						Kout(i3 * 3 + x, i3 * 3 + y) += K[2](x, y);
					}
				}
			}
		}
	}

	void FemPhysicsMatrixFree::SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, real pressure)
	{
		mAppliedPressure = pressure;
		mTractionSurface.resize(triangleList.size());

		for (uint32 i = 0; i < triangleList.size(); i++)
		{
			mTractionSurface[i] = triangleList[i];
		}
	}

} // namespace FEM_SYSTEM