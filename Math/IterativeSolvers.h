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

#ifndef ITERATIVE_SOLVERS_H
#define ITERATIVE_SOLVERS_H

// References:
// [Shewchuck] Shewchuk, J.R., "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain"

// some utilities for treating std::vector as a linear algebra vector
template <typename Real>
std::vector<Real> operator -(const std::vector<Real>& a, const std::vector<Real>& b)
{
	size_t n = a.size();
	std::vector<Real> r(n);
	for (size_t i = 0; i < n; i++)
		r[i] = a[i] - b[i];
	return r;
}

template <typename Real>
std::vector<Real> operator +(const std::vector<Real>& a, const std::vector<Real>& b)
{
	size_t n = a.size();
	std::vector<Real> r(n);
	for (size_t i = 0; i < n; i++)
		r[i] = a[i] + b[i];
	return r;
}

template <typename Elem>
std::vector<Elem> operator *(float s, const std::vector<Elem>& v)
{
	size_t n = v.size();
	std::vector<Elem> r(n);
	for (size_t i = 0; i < n; i++)
		r[i] = s * v[i];
	return r;
}

template <typename Elem>
std::vector<Elem> operator *(double s, const std::vector<Elem>& v)
{
	size_t n = v.size();
	std::vector<Elem> r(n);
	for (size_t i = 0; i < n; i++)
		r[i] = s * v[i];
	return r;
}

// inner product
template <typename Elem>
float operator *(const std::vector<Elem>& a, const std::vector<Elem>& b)
{
	ASSERT(a.size() == b.size());
	float sum = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		float dot = a[i] * b[i];
		sum += dot;
	}
	return sum;
}

template <typename TE, typename TR>
TR InnerProduct(const std::vector<TE>& a, const std::vector<TE>& b)
{
	ASSERT(a.size() == b.size());
	TR sum = 0;
	for (size_t i = 0; i < a.size(); i++)
		sum += a[i] * b[i];
	return sum;
}

template <typename T1, typename T2, typename TR>
TR InnerProduct(const std::vector<T1>& a, const std::vector<T2>& b)
{
	ASSERT(a.size() == b.size());
	TR sum = 0;
	for (size_t i = 0; i < a.size(); i++)
		sum += a[i] * b[i];
	return sum;
}

// Matrix free methods
template <typename Real, typename Vector, class MatrixMultiplier>
int SolveConjugateGradientMF(const MatrixMultiplier& A, const std::vector<Vector>& b, std::vector<Vector>& x, int numIters, float eps = 1e-3f, bool usePrecond = false)
{
	// algorithm taken from [Shewchuk]

	std::vector<Vector> r = b; // allocation!
	// compute the full residual only if x != 0
	std::vector<Vector> df(b.size()); // allocation!
	if (InnerProduct<Vector, Real>(x, x) != 0)
	{
		A.MatrixVectorMultiply(x, df);
		r = b - df;
	}
	std::vector<Vector> d = r; // allocation!

	// initialize vector used for preconditioning
	std::vector<Vector> pcv;
	std::vector<Real> diagInv;
	if (usePrecond)
	{
		// allocate only if needed
		pcv.resize(b.size());
		diagInv.resize(b.size());

		// pre-compute inverse diagonal values
		for (size_t i = 0; i < diagInv.size(); i++)
		{
			// prepare unit vector
			std::fill(pcv.begin(), pcv.end(), Vector(0));
			pcv[i] = Vector(1);
			// compute A-dot product
			A.MatrixVectorMultiply(pcv, df);
			Real val = InnerProduct<Vector, Real>(pcv, df);
			diagInv[i] = Real(1) / val;

			d[i] = r[i] * diagInv[i];
		}
	}
	std::vector<Vector>& s = usePrecond ? pcv : r; // just a shell for r if not preconditioned

	Real delta = InnerProduct<Vector, Real>(r, d);
	if (delta == 0)
		return 0;
	Real delta0 = delta;

	for (int iter = 0; iter < numIters; iter++)
	{
		A.MatrixVectorMultiply(d, df);
		Real alpha = delta / InnerProduct<Vector, Real>(d, df);
		x = x + alpha * d;
		// TODO: move the calculations below to the top of the loop
#ifdef EXACT_RESIDUAL
		if (iter % 50 == 0)
		{
			A.MatrixVectorMultiply(x, df);
			r = b - df;
		}
		else
#endif
		r = r - alpha * df;

		if (usePrecond)
		{
			for (size_t i = 0; i < diagInv.size(); i++)
				s[i] = r[i] * diagInv[i];
		}

		Real delta1 = InnerProduct<Vector, Real>(r, s);
		Real beta = delta1 / delta;
		delta = delta1;

		// the convergence criterion from [Shewchuk]
		if (abs(delta) <= eps * eps * abs(delta0))
		{
			// convergence criterion met
			return iter + 1;
		}
		d = s + beta * d;
	}
	return numIters;
}

template <typename Real, typename Vector, class MatrixMultiplier>
int SolveConjugateResidualMF(const MatrixMultiplier& A, const std::vector<Vector>& b, std::vector<Vector>& x, int numIters, float eps = 1e-3f, bool usePrecond = false)
{
	// algorithm taken from [Shewchuk]

	std::vector<Vector> r = b; // allocation!
	// compute the full residual only if x != 0
	std::vector<Vector> q(b.size()); // allocation!
	if (InnerProduct<Vector, Real>(x, x) != 0)
	{
		A.MatrixVectorMultiply(x, q); // q = Ax
		r = b - q;
	}

	// initialize vector used for preconditioning
	std::vector<Vector> pcv;
	std::vector<Real> diagInv;
	if (usePrecond)
	{
		// allocate only if needed
		pcv.resize(b.size());
		diagInv.resize(b.size());

		// pre-compute inverse diagonal values
		for (size_t i = 0; i < diagInv.size(); i++)
		{
			// prepare unit vector
			std::fill(pcv.begin(), pcv.end(), Vector(0));
			pcv[i] = Vector(1);
			// compute A-dot product
			A.MatrixVectorMultiply(pcv, q);
			Real val = InnerProduct<Vector, Real>(pcv, q);
			diagInv[i] = Real(1) / val;

			r[i] = r[i] * diagInv[i];
		}
	}

	std::vector<Vector> d = r; // allocation!	
	A.MatrixVectorMultiply(r, q); // q = Ar

	std::vector<Vector> p(b.size()); // allocation!
	A.MatrixVectorMultiply(d, p); // p = Ad
	std::vector<Vector>& s = usePrecond ? pcv : p; // just a shell for p if not preconditioned

	Real delta = InnerProduct<Vector, Real>(r, q);
	if (delta == 0)
		return 0;
	Real delta0 = InnerProduct<Vector, Real>(r, r);

	for (int iter = 0; iter < numIters; iter++)
	{
		if (usePrecond)
		{
			for (size_t i = 0; i < diagInv.size(); i++)
				s[i] = p[i] * diagInv[i];
		}

		Real alpha = delta / InnerProduct<Vector, Real>(p, s);
		x = x + alpha * d;
		r = r - alpha * s;

		A.MatrixVectorMultiply(r, q); // q = Ar

		Real delta1 = InnerProduct<Vector, Real>(r, q);
		Real beta = delta1 / delta;
		delta = delta1;

		// the convergence criterion from [Shewchuk]
		if (abs(InnerProduct<Vector, Real>(r, r)) <= eps * eps * abs(delta0))
		{
			// convergence criterion met
			return iter + 1;
		}
		
		d = r + beta * d;
		p = q + beta * p;
	}
	return numIters;
}

#endif // ITERATIVE_SOLVERS_H