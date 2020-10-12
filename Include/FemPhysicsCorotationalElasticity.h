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

#ifndef FEM_PHYSICS_ANY_ORDER_COROTATIONAL_ELASTICITY_H
#define FEM_PHYSICS_ANY_ORDER_COROTATIONAL_ELASTICITY_H

#include "FemPhysicsLinear.h"
#include "TetrahedralMesh.h"
#include "LinearSolver.h"
#include <memory>

namespace FEM_SYSTEM
{
	class FemPhysicsCorotationalElasticity : public FemPhysicsLinear
	{
	public:
		FemPhysicsCorotationalElasticity(std::vector<Tet>& tetrahedra,
			std::vector<Node>& allNodes, const FemConfig& config);

		void Step(real dt) override;
		void SolveEquilibrium(float);

	private:
		void StepExplicit(real h);
		void StepImplicit(real h);
		void AssembleStiffnessMatrix();
		void CacheLocalStiffnessMatrices();
		void ComputeRotations();
		void ComputeElasticForces(EigenVector& elasticForce) const;

	private:
		std::vector<EigenMatrix> mCachedLocalStiffnessMatrix;
		SparseMatrix mStiffnessMatrix; // stiffness matrix
		EigenMatrix mDenseK; // dense stiffness matrix
		std::vector<Matrix3R> mCachedRotationMatrix;
		EigenVector mElasticForce;
		LinearSolver mExplicitMassSolver;
	};


} // namespace FEM_SYSTEM

#endif // FEM_PHYSICS_ANY_ORDER_COROTATIONAL_ELASTICITY_H
