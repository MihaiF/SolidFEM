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

#ifndef FEM_PHYSICS_MIXED
#define FEM_PHYSICS_MIXED

#include "FemPhysicsLinearIncompressible.h"

namespace FEM_SYSTEM
{
	class PrimalProblem;

	class FemPhysicsMixed : public FemPhysicsLinearIncompressible
	{
	public:
		FemPhysicsMixed(std::vector<Tet>& tetrahedra,
			std::vector<Node>& allNodes, const FemConfig& config);

		void Step(real dt) override;

	private:
		void SolveSaddlePoint();
		void SolveSaddlePointLS();

		real GetCurrentVolumeErrors(EigenVector& errors, bool verbose = false);
		void AssembleGeometricStiffnessMatrix(const EigenVector& p, SparseMatrix& Kgeom) const;
		void AssembleGeometricStiffnessMatrixFD(const EigenVector& p, EigenMatrix& Kgeom);
		void ComputeGradients(std::vector<Vector3R>& r, const EigenVector& dx, int material = -1);

		template<class MATRIX>
		void ComputeStiffnessMatrix(MATRIX& K, bool update = true);
		EigenVector ComputeStdRhs(int material = -1);
		EigenVector ComputePosRhs();

		void AssembleJacobianMatrixFD(EigenMatrix& J);

		uint32 GetNumContacts() const
		{
			return 0;
		}

	private:
		real mTimeStep;
		EigenVector mContactMultipliers, mContactDepth;
		bool mCondenseContact = false;
		real mStepLength;
		EigenVector mVolStretches;

		friend class SaddlePointProblem;
	};

} // namespace FEM_SYSTEM

#endif // FEM_PHYSICS_MIXED