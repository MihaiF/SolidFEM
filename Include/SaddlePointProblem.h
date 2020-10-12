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

#pragma once

#include "MixedProblems.h"

namespace FEM_SYSTEM
{
	class SaddlePointProblem : public BaseProblem
	{
	public:
		enum SolverType
		{
			SPP_ST_DIRECT,
			SPP_ST_ITERATIVE,
			SPP_ST_SCHUR_COMPLIANCE,
			SPP_ST_SCHUR_MASS,
		};
	public:
		SaddlePointProblem(FemPhysicsMixed& source, int solverType) : BaseProblem(source), mSolverType(solverType)
		{
			if (mSolverType == SPP_ST_SCHUR_COMPLIANCE)
			{
				EigenMatrix C(mFPM.mVolComplianceMatrix);
				mCinv = C.inverse().sparseView();
			}
		}
		EigenVector ComputeRhs(const EigenVector& solution);
		void ComputeSystemMatrix(const EigenVector& s, const EigenVector& y, EigenMatrix& K);
		void ComputeSystemMatrix(const EigenVector& s, const EigenVector& y, SparseMatrix& K);
		EigenVector SolveLinearSystem(EigenMatrix& K, const EigenVector& rhs, const EigenVector& s, const EigenVector& y);
		EigenVector SolveLinearSystem(SparseMatrix& K, const EigenVector& rhs, const EigenVector& s, const EigenVector& y);
		real MeritResidual(const EigenVector& rhs) const;

	private:
		SparseMatrix mCinv; // inverse compliance matrix
		int mSolverType = SPP_ST_DIRECT;
		EigenMatrix mKinv; // inverse stiffness matrix
	};
}
