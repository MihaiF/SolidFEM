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

#ifndef FEM_PHYSICS_H
#define FEM_PHYSICS_H
#include "FemDataStructures.h"
#include "FemPhysicsBase.h"
#include <Math/IterativeSolvers.h> // for InnerProduct

namespace FEM_SYSTEM
{
	class FemPhysicsMatrixFree;

	class FemPhysicsMatrixFree : public FemPhysicsBase
	{
		public:
			enum NodalForcesType
			{
				NFT_DIRECT,
				NFT_FINITE_VOLUME,
			};

			struct Config
			{
				NodalForcesType mNodalForcesComputation = NFT_DIRECT; // not used
				NonlinearSolverType mSolver = NST_NEWTON;
				float mDescentRate = 1e-5f;
				bool mOptimizer = false;
			};

		public:
			FemPhysicsMatrixFree(std::vector<Tet>& tetrahedra,
				std::vector<Node>& allNodes, const FemConfig& config);

			void Step(real dt) override;
			
			const std::vector<Node>& GetNodes() const;
			std::vector<Node>& GetNodes(); // should we allow this?
			const std::vector<Tetrahedron>& GetTetrahedra() const;
			void UpdatePositions(std::vector<Node>& nodes) override;
			Vector3R GetDeformedPosition(uint32 idx, bool shuffle = true) const override { return nodes[shuffle ? mReshuffleMap[idx] : idx].pos; }
			Vector3R GetInitialPosition(uint32 idx) const override { return nodes[mReshuffleMap[idx]].pos0; } // Careful! These get over-written by implicit integrator
			void SetDeformedPosition(uint32 idx, Vector3R val, bool shuffle = true) override { nodes[shuffle ? mReshuffleMap[idx] : idx].pos = val; }
			Vector3R GetVelocity(uint32 idx) const override { return nodes[mReshuffleMap[idx]].vel; }
			uint32 GetNumNodes() const override { return (uint32)nodes.size(); }
			uint32 GetNumLocalNodes() const override { return 4; }
			uint32 GetNumFreeNodes() const override { ASSERT(mNumBCs < GetNumNodes()); return GetNumNodes() - mNumBCs; }
			uint32 GetNumElements() const override { return (uint32)tets.size(); }
			bool IsNodeFixed(uint32 i) const override { ASSERT(i < GetNumNodes()); return nodes[mReshuffleMap[i]].invMass == 0; }
			Matrix3R GetBarycentricJacobianMatrix(uint32 e) const override { return tets[e].Xtr; }
			real GetElementInitialVolume(uint32 e) const override { return tets[e].vol; }
			uint32 GetGlobalIndex(uint32 e, uint32 l) const override { return tets[e].i[l]; }
			uint32 GetGlobalOriginalIndex(uint32 e, uint32 l) const override { return originalTets[e].i[l]; }
			real GetTotalVolume() const;

			// Krylov solvers helpers
			void MatrixVectorMultiply(const std::vector<Vector3R>& in, std::vector<Vector3R>& out) const; // used for matrix free CG
			void ComputeGradients(std::vector<Vector3R>& r);
			void UpdatePosAndComputeGradients(const std::vector<Vector3R>& pos, std::vector<Vector3R>& r);
			real DotProduct(const std::vector<Vector3R>& a, const std::vector<Vector3R>& b) const { return InnerProduct<Vector3R, real>(a, b); }

			real ComputeEnergy(); // basically the minimization objective
			real ComputeEnergy(int level);

			// Newton solver helpers
			EigenVector ComputeRhs(const EigenVector& solution);			
			template<class MATRIX> void ComputeSystemMatrix(const EigenVector& s, const EigenVector& y, MATRIX& K);
			EigenVector SolveLinearSystem(EigenMatrix& K, const EigenVector& rhs, const EigenVector& s, const EigenVector& y);
			EigenVector SolveLinearSystem(SparseMatrix& K, const EigenVector& rhs, const EigenVector& s, const EigenVector& y);
			real MeritResidual(const EigenVector& rhs);

			Vector3R GetTotalDisplacement(uint32 i) const { return (nodes[i + mNumBCs].pos - nodes[i + mNumBCs].pos0); }
			void SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, real pressure) override;

			real ComputeSpringEnergy(Cable& cable);
			real ComputeSpringEnergy();

		private:
			// traditional FEM
			void SubStep(real h);
			bool Solve();

			// nonlinear static solvers
			void SolveNewtonCG();
			bool SolveNewton();
			bool SolveNewtonLS();
			void SolveGradientDescent(real alpha);
			void SolveNonlinearSteepestDescent();
			void SolveNonlinearConjugateGradient();

			void ComputeDeformationGradient(uint32 e, Matrix3R& F) const override;

			void ComputePressureForces(Vector3Array& fout, EigenMatrix& K) const;

			void ReshuffleFixedNodes();

			void AddCable(const std::vector<SpringNode>& cable, const Vector3Array& pos, real restLength, real stiffness, real damping, real actuation) override;

		protected:
			void ComputeForceDifferential(const std::vector<Vector3R>& dx, std::vector<Vector3R>& df) const;
			void BuildMassMatrix();

		protected:
			std::vector<Node> nodes;
			std::vector<Tetrahedron> tets, originalTets;
			real ed, nud;
			Matrix3R Ed;
			bool hasCollisions;
			std::vector<Vector3R> disp; // displacements vector a.k.a u
			std::vector<Vector6, Eigen::aligned_allocator<Vector6>> lambdaAcc; // accumulated lambdas

			Config mConfig;

			real mTimeStep = 0;

			std::vector<Vector3R> mForces;
			SparseMatrix mMassMatrix;

			mutable std::vector<Matrix3R> dH;

			uint32 mNumSpringNodes = 0;

			friend class FemSystem;
			friend class ElasticProblem;
	};

	inline const std::vector<Node>& FemPhysicsMatrixFree::GetNodes() const
	{
		return nodes;
	}

	inline std::vector<Node>& FemPhysicsMatrixFree::GetNodes()
	{
		return nodes;
	}

	inline const std::vector<Tetrahedron>& FemPhysicsMatrixFree::GetTetrahedra() const
	{
		return tets;
	}

}
#endif // FEM_PHYSICS_MATRIX_FREE_H