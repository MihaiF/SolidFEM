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

#ifndef NONLINEAR_CG_H
#define NONLINEAR_CG_H

namespace FEM_SYSTEM
{
	template<class PROBLEM, class VECTOR>
	class NonlinearConjugateGradientMinimizer
	{
	public:
		void Solve(PROBLEM& problem, uint32 size, VECTOR& solution);

	public:
		uint32 mOuterIterations = 200;
		uint32 mInnerIterations = 10;
		real mAbsResidualThreshold = 1e-4;
		bool mVerbose = false;
	};
	
	template<class PROBLEM, class VECTOR>
	void NonlinearConjugateGradientMinimizer<PROBLEM, VECTOR>::Solve(PROBLEM& problem, uint32 size, VECTOR& solution)
	{
		VECTOR r(size); // allocation!
		VECTOR df(size); // allocation!
		VECTOR d(size); // allocation!

		const float eps = 1e-2f;
		const float tol = 1e-2f;

		problem.ComputeGradients(r);

		d = r;
		real delta = problem.DotProduct(r, r);
		if (mVerbose)
			Printf("%g\n", sqrt(delta));
		real delta0 = delta;
		if (delta == 0)
			return;

		// nonlinear CG
		for (uint32 iter = 0; iter < mOuterIterations; iter++)
		{
			// line search step (Newton-Raphson)
			real deltaD = problem.DotProduct(d, d);
			for (uint32 innerIter = 0; innerIter < mInnerIterations; innerIter++)
			{
				problem.MatrixVectorMultiply(d, df);
				real alpha = problem.DotProduct(r, d) / problem.DotProduct(d, df);

				// update positions
				for (size_t i = 0; i < size; i++)
				{
					solution[i] += alpha * d[i];
				}

				problem.UpdatePosAndComputeGradients(solution, r);

				// line search convergence criterion
				if (alpha * alpha * deltaD < tol * tol)
				{
					//Printf("Inner convergence reached after %d iterations.\n", innerIter + 1);
					break;
				}
			}
			real delta1 = problem.DotProduct(r, r);
			real beta = delta1 / delta;
			delta = delta1;
			if (mVerbose)
				Printf("%g\n", sqrt(delta));
			if (delta < mAbsResidualThreshold * mAbsResidualThreshold)
			{
				if (mVerbose)
					Printf("Absolute convergence reached after %d iterations.\n", iter + 1);
				break;
			}
			if (delta <= eps * eps * delta0)
			{
				//Printf("Convergence reached after %d iterations.\n", iter + 1);
				break;
			}
			d = r + beta * d;
		}
	}

}

#endif // NONLINEAR_CG_H