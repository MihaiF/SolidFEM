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

#ifndef FEM_PHYSICS_LINEAR_EXPLICIT_H
#define FEM_PHYSICS_LINEAR_EXPLICIT_H

#include "FemPhysicsLinear.h"
#include "LinearSolver.h"

namespace FEM_SYSTEM
{
	class FemPhysicsLinearElasticity : public FemPhysicsLinear
	{
	public:
		FemPhysicsLinearElasticity(std::vector<Tetrahedron>& tetrahedra,
			std::vector<Node>& nodes, const FemConfig& config);
		void Step(real dt) override;
		void SolveEquilibrium(float) override;

	private:
		void StepExplicit(real h);
		void StepImplicit(real h);
		void AssembleStiffnessMatrix();

	protected:
		// linear elasticity helpers
		void ComputeLocalStiffnessMatrix(uint32 i, EigenMatrix& Klocal);
		void ComputeLocalStiffnessMatrixBB1(uint32 i, EigenMatrix& Klocal);
		void ComputeDyadicMatrixBB2(uint32 i, uint32 j, Vector3R y[4], Matrix3R& L);
		void ComputeLocalStiffnessMatrixBB2(uint32 i, EigenMatrix& Klocal);

	protected:
		EigenMatrix mStiffnessMatrix; // stiffness matrix
		Eigen::SparseLU<SparseMatrix> mSparseLU; // the LU decomposition of M
		
		// for implicit integration
		bool mHaveSysMatrix = false;
		LinearSolver mSolver;

		friend class FemTester;
	};

} // namespace FEM_SYSTEM

#endif // !FEM_PHYSICS_LINEAR_EXPLICIT_H
