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

#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include "LinearSolver.h"
#include <Engine/Types.h>
#include <Engine/Utils.h>

#include <iostream>

//namespace FEM_SYSTEM
//{
	enum VerboseLevel
	{
		VL_NONE = 0,
		VL_RESIDUAL,
		VL_MINIMUM,
		VL_DETAILED,
		VL_DEBUG,
	};

	enum LSCondition
	{
		LSC_ARMIJO,
		LSC_ENGINEERING,
	};

	// A simple Newton solver with no line search and no convergence criterion
	template<class PROBLEM, class MATRIX = EigenMatrix>
	class NewtonSolver
	{
	public:
		void Solve(PROBLEM& problem, uint32 size, EigenVector& solution);
		void ComputeJacobianFD(PROBLEM& problem, uint32 size, EigenVector& solution, MATRIX& K);

	public:
		int mNumIterations = 200;
		real mResidualThreshold = 1;
		bool mDebug = true;
		int mVerbose = 1;
		LinearSolverType mSolverType = LST_LU_PARDISO;
		bool mUseFiniteDiff = false;
	};

	// A naive back-tracking Newton solver of my own (halving the step)
	template<class PROBLEM, class MATRIX = EigenMatrix>
	class NewtonSolverBackTrack : public NewtonSolver<PROBLEM, MATRIX>
	{
	public:
		bool Solve(PROBLEM& problem, uint32 size, EigenVector& solution);

	public:
		int mLSCondition = LSC_ENGINEERING;
		real mAlpha = 1e-4;
		real mBackTrackFactor = 0.5;
		int mMaxBisection = 10;
		bool mUseProblemSolver = false;

		bool mUseBFGS = false;
		MATRIX mB; // B matrix in BFGS
		EigenMatrix mH; // H matrix in BFGS
	};

	template<class PROBLEM, class MATRIX>
	void NewtonSolver<PROBLEM, MATRIX>::ComputeJacobianFD(PROBLEM& problem, uint32 size, EigenVector& solution, MATRIX& K)
	{
		real eps = 1e-5; 
		K.resize(size, size);
		EigenVector rhs0 = problem.ComputeRhs(solution);
		for (uint32 j = 0; j < size; j++)
		{
			solution(j) += eps;
			EigenVector rhs = problem.ComputeRhs(solution);
			//K.block(0, j, size, 1) = -(rhs - rhs0) / eps;
			auto d = -(rhs - rhs0) / eps;
			for (uint32 i = 0; i < size; i++)
				K.coeffRef(i, j) = d(i);
			solution(j) -= eps;
		}
	}

	template<class PROBLEM, class MATRIX>
	void NewtonSolver<PROBLEM, MATRIX>::Solve(PROBLEM& problem, uint32 size, EigenVector& solution)
	{
		// working vars
		MATRIX K;
		EigenVector s(0), y(0); // dummy vars

		// assemble rhs
		EigenVector rhs = problem.ComputeRhs(solution);
		real res = rhs.norm();
		if (mVerbose)
		{
			Printf("%g\n", res);
		}

		// Newton steps (outer loop)
		bool converged = false;
		for (int iter = 0; iter < mNumIterations; iter++)
		{
			// assemble system matrix K
			if (mUseFiniteDiff)
				ComputeJacobianFD(problem, size, solution, K);
			else
				problem.ComputeSystemMatrix(s, y, K); // we assume ComputeRhs has already been called, so we don't need to send the current solution once more			
			if (mVerbose > VL_DETAILED)
				std::cout << "K" << std::endl << K << std::endl;

			// solve K * dx = b using Eigen
			LinearSolver solver;
			solver.Init(K, mSolverType);
			EigenVector delta = solver.Solve(rhs);

			solution = solution + delta;

			// assemble rhs
			rhs = problem.ComputeRhs(solution); // this also sends the final solution to the caller problem

			// compute the residual
			real res = rhs.norm();
			if (mVerbose)
			{
				Printf("%g\n", res);
			}
			if (mDebug && isnan(res))
			{
				if (mVerbose)
					Printf("nan rhs\n");
				break;
			}
			if (mDebug && problem.CheckForInversion())
			{
				if (mVerbose)
					Printf("Inverted elements\n");
			}
			if (res < mResidualThreshold)
			{
				if (mVerbose)
					Printf("Absolute residual convergence at iteration %d\n", iter + 1);
				converged = true;
				break;
			}
		}

		if (mVerbose && !converged)
			Printf("The Newton solver has not converged within the max number of iterations. Consider reducing the time step.\n");
	}

	template<class PROBLEM, class MATRIX>
	bool NewtonSolverBackTrack<PROBLEM, MATRIX>::Solve(PROBLEM& problem, uint32 size, EigenVector& solution)
	{
		// working vars
		EigenVector delta, solutionOld, rhs, newRhs;
		EigenVector s(0), y(0); // used by BFGS

		// assemble rhs
		rhs = problem.ComputeRhs(solution);
		real res = rhs.norm(); // residual
		if (this->mVerbose)
			Printf("%g\n", res);

		// Newton steps (outer loop)
		bool converged = false;
		for (int iter = 0; iter < this->mNumIterations; iter++)
		{

			// compute BFGS approximation
			if (iter > 0)
			{
				// compute the BFGS vectors anyway, in case the problem wants to use them
				s = solution - solutionOld; // or damping * delta
				y = newRhs - rhs;
				if (mUseBFGS)
				{
					mB = mB - (mB * s * s.transpose() * mB) / (s.transpose() * mB * s) + (y * y.transpose()) / (y.transpose() * s);
					//mH = mH - (mH * y * y.transpose() * mH) / (y.transpose() * mH * y) + (s * s.transpose()) / (y.transpose() * s);
				}
			}

			if (!mUseBFGS || iter == 0)
			{
				// assemble system matrix K
				problem.ComputeSystemMatrix(s, y, mB); // we assume ComputeRhs has already been called, so we don't need to send the current solution once more
			}

			if (!mUseProblemSolver)
			{
				// solve K * dx = b using Eigen
				LinearSolver solver;
				solver.Init(mB, this->mSolverType);
				delta = solver.Solve(rhs);
			}
			else
				delta = problem.SolveLinearSystem(mB, rhs, s, y);

			// save old solution
			solutionOld = solution;

			// try to find a suitable step length - "line search"
			real damping = 1;
			real Rold = -delta.transpose() * rhs;
			real fOld = problem.MeritResidual(rhs);
			auto grad = -rhs; // when -rhs is the gradient of the objective
			real slope = delta.transpose() * grad;
			if (slope > 0 && this->mVerbose > VL_MINIMUM)
				Printf("positive slope %e\n", slope); // this is mostly for when minimizing the norm
			real minStep = 0;
			real maxStep = 1;
			int bisectionIterations = 0;
			bool giveup = false;
			while (true)
			{
				// update the solution
				solution = solutionOld + damping * delta;

				// compute new RHS (minus function)
				newRhs = problem.ComputeRhs(solution);

				// check if the solution is acceptable
				bool cond = true;
				if (mLSCondition == LSC_ARMIJO)
				{
					real fNew = problem.MeritResidual(newRhs);
					cond = fNew < fOld + mAlpha * damping * slope; // Armijo condition
				}
				else
				{
					real Rnew = -delta.transpose() * newRhs;
					cond = Rnew < 0.5 * Rold; // Bonet-Wood condition
				}
				if (!cond) // NR condition
				{
					maxStep = damping;
					bisectionIterations++;						
				}
				else
				{
					if (bisectionIterations > 10 // reached the bisection limit
						|| damping == 1) // or full step works
					{
						if (this->mVerbose > VL_MINIMUM)
							Printf("acceptable solution found\n");
						break;
					}
					else
					{
						minStep = damping;
						bisectionIterations++;
					}
				}

				damping = 0.5 * (minStep + maxStep); // TODO: use largest value from the roots
				if (this->mVerbose > VL_MINIMUM)
					Printf("step length %g\n", damping);
				if (damping < 1e-10)
				{
					if (this->mVerbose > VL_MINIMUM)
						Printf("step reduced too much - quitting\n");
					giveup = true;
					break;
				}
			} // end of back-track iteration
			if (this->mVerbose > VL_RESIDUAL)
			{
				if (damping < 1)
					Printf("step length: %.2g\n", damping);
			}

			// compute the new residual
			real newRes = newRhs.norm();
			// compute the relative change
			real rel = (res - newRes) / res;
			if (this->mVerbose)
			{
				if (this->mVerbose > VL_MINIMUM)
					Printf("residual %d: %g / %g (%e)\n", iter, newRes, res, rel);
				else
					Printf("%g\n", newRes);
			}

			// quit conditions
			if (this->mDebug && isnan(newRes))
			{
				if (this->mVerbose)
					Printf("nan rhs\n");
				solution = solutionOld;
				break;
			}
			if (giveup)
			{
				solution = solutionOld;
				break;
			}
			if (newRes < this->mResidualThreshold)
			{
				if (this->mVerbose > VL_RESIDUAL)
					Printf("Absolute residual convergence at iteration %d\n", iter + 1);
				converged = true;
				break;
			}
			if (abs(rel) < 1e-10)
			{
				if (this->mVerbose)
					Printf("The residual is decreasing too slow. Quitting at iteration %d.\n", iter + 1);
				break;
			}

			// update the RHS for the next step
			rhs = newRhs;
			res = newRes;
		} // end of Newton iteration

		if (this->mVerbose && !converged)
			Printf("The Newton solver has not converged within the max number of iterations. Consider reducing the time step.\n");

		return converged;
	}

//} // namespace FEM_SYSTEM

#endif // NEWTON_SOLVER_H
