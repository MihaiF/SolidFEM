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

#include "PolarDecomposition.h"

namespace FEM_SYSTEM
{
	void ComputePolarDecomposition(const Matrix3R& F, Matrix3R& R, Matrix3R& U)
	{
		PolarDecomposition::QRPolarDecomposition(F, R);
		Matrix3R V = !R * F;
		U = 0.5f * (V + !V); // make U symmetric from V

		//if (!PolarDecomposition::EigenPolarDecomposition(F, R, U) || R.Determinant() < 0)
		//{
		//	Printf("Polar decomposition failed\n");
		//}
	}

	// Function below adapted from OpenTissue: https://github.com/erleben/OpenTissue/blob/master/OpenTissue/core/math/math_eigen_system_decomposition.h
	/**
		* Eigen System Decomposition.
		*
		* @param A        Matrix to find eigenvectors and eigenvalues of.
		* @param V        Upon return the columns of this matrix contains the
		*                 eigenvectors. IMPORTANT: V may not form a right-handed
		*                 coordinate system (ie. V might not be a rotation matrix). If
		*                 a rotation matrix is needed one should compute the determinant
		*                 of V and flip the sign of one of V's columns in case the
		*                 determinant is negative. That is:
		*
		*                    if(det(V)<0) {  V(0,0) =- V(0,0); V(1,0) =- V(1,0); V(2,0) =- V(2,0); }
		*
		*
		* @param diag     Upon return this vector contains the
		*                 eigenvalues, such that entry 0 correpsonds to
		*                 eigenvector 0 and so on.
		*/
	void EigenDecomposition(Matrix3R const & A, Matrix3R & V, Vector3R & diag)
	{
		Vector3R sub_diag;

		sub_diag.SetZero();
		diag.SetZero();

		V = A;

		real const & fM00 = V(0, 0);
		real fM01 = V(0, 1);
		real fM02 = V(0, 2);
		real const & fM11 = V(1, 1);
		real const & fM12 = V(1, 2);
		real const & fM22 = V(2, 2);

		diag[0] = fM00;
		sub_diag[2] = 0;
		if (fM02 != 0)
		{
			real fLength = sqrt(fM01*fM01 + fM02 * fM02);
			real fInvLength = 1 / fLength;
			fM01 *= fInvLength;
			fM02 *= fInvLength;
			real fQ = 2*fM01*fM12 + fM02 * (fM22 - fM11);
			diag[1] = fM11 + fM02 * fQ;
			diag[2] = fM22 - fM02 * fQ;
			sub_diag[0] = fLength;
			sub_diag[1] = fM12 - fM01 * fQ;
			V(0, 0) = 1;
			V(0, 1) = 0;
			V(0, 2) = 0;
			V(1, 0) = 0;
			V(1, 1) = fM01;
			V(1, 2) = fM02;
			V(2, 0) = 0;
			V(2, 1) = fM02;
			V(2, 2) = -fM01;
		}
		else
		{
			diag[1] = fM11;
			diag[2] = fM22;
			sub_diag[0] = fM01;
			sub_diag[1] = fM12;
			V(0, 0) = 1;
			V(0, 1) = 0;
			V(0, 2) = 0;
			V(1, 0) = 0;
			V(1, 1) = 1;
			V(1, 2) = 0;
			V(2, 0) = 0;
			V(2, 1) = 0;
			V(2, 2) = 1;
		}

		const int max_iterations = 32;
		const int dim = 3;
		for (int i0 = 0; i0 < dim; ++i0)
		{
			int i1;
			for (i1 = 0; i1 < max_iterations; ++i1)
			{
				int i2;
				for (i2 = i0; i2 <= dim - 2; ++i2)
				{
					real fTmp = fabs(diag(i2)) + fabs(diag(i2 + 1));
					if (fabs(sub_diag(i2)) + fTmp == fTmp)
						break;
				}
				if (i2 == i0)
					break;
				real fG = (diag(i0 + 1) - diag(i0)) / (2*  sub_diag(i0));
				real fR = sqrt(fG*fG + 1);
				if (fG < 0)
					fG = diag(i2) - diag(i0) + sub_diag(i0) / (fG - fR);
				else
					fG = diag(i2) - diag(i0) + sub_diag(i0) / (fG + fR);

				real fSin = 1;
				real fCos = 1;
				real fP = 0;

				for (int i3 = i2 - 1; i3 >= i0; --i3)
				{
					real fF = fSin * sub_diag(i3);
					real fB = fCos * sub_diag(i3);
					if (fabs(fF) >= fabs(fG))
					{
						fCos = fG / fF;
						fR = sqrt(fCos*fCos + 1);
						sub_diag(i3 + 1) = fF * fR;
						fSin = 1 / fR;
						fCos *= fSin;
					}
					else
					{
						fSin = fF / fG;
						fR = sqrt(fSin*fSin + 1);
						sub_diag(i3 + 1) = fG * fR;
						fCos = 1 / fR;
						fSin *= fCos;
					}
					fG = diag(i3 + 1) - fP;
					fR = (diag(i3) - fG)*fSin + 2*fB*fCos;
					fP = fSin * fR;
					diag(i3 + 1) = fG + fP;
					fG = fCos * fR - fB;
					for (int i4 = 0; i4 < dim; ++i4)
					{
						fF = V(i4, i3 + 1);
						V(i4, i3 + 1) = fSin * V(i4, i3) + fCos * fF;
						V(i4, i3) = fCos * V(i4, i3) - fSin * fF;
					}
				}
				diag(i0) -= fP;
				sub_diag(i0) = fG;
				sub_diag(i2) = 1;
			}
			if (i1 == max_iterations)
				break;
		}
	}

	// Function below adapted from OpenTissue: https://github.com/erleben/OpenTissue/blob/master/OpenTissue/core/math/math_polar_decomposition.h
	//* Polar Decomposition of matrix A (as described by Etzmuss et. al in ``A Fast Finite Solution for Cloth Modelling'')
	//*
	//*   A = R S
	//*
	//* Where R is a special orthogonal matrix (R*R^T = I and det(R)=1)
	//*
	//*   S^2 = A^T A
	//*
	//* let d be vector of eigenvalues and let v_0,v_1, and v_2 be corresponding eigenvectors of (A^T A), then
	//*
	//*   S = sqrt(d_0) v_0 * v_0^T + ... + sqrt(d_2) v_2 * v_2^T
	//*
	//* Now compute
	//*
	//*  R = A * S^-1
	bool PolarDecomposition::EigenPolarDecomposition(Matrix3R const & A, Matrix3R & R, Matrix3R & S)
	{
		// A = U D V^T       // A is square so: thin SVD = SVD
		// A = (U V^T)  (V D V^T)
		//   =    R        S
		//
		//Notice S is symmetric  R should be orthonormal?
		//
		//  proof of
		//      S =  sqrt( A^T A )
		//
		//      start by
		//
		//            S * S = A^T A
		// V D V^T V D V^T  =   V D U^T  U D V^T
		//       V D D V^T  =   V D D V^T
		//Assume
		//A = R S
		//pre-multiply and use assumption then
		//  A^T A = A^T R S
		//  A^T A = S^T S = S S
		//  last step used S is symmetric

		// Find the eigenvalues of A^T A using Eigen library
		Eigen::Matrix<real, 3, 3> Aeig = Matrix3ToEigen<Eigen::Matrix<real, 3, 3>>(A);
		auto S2 = Aeig.transpose() * Aeig;

		Eigen::EigenSolver< Eigen::Matrix<real, 3, 3>> eigenSolver(S2);
		auto V = eigenSolver.eigenvectors().real();
		auto d = eigenSolver.eigenvalues().real();
		
		//Matrix3R V;
		//Vector3R d;
		//EigenDecomposition(!A * A, V, d);

		//--- Test if all eigenvalues are positive (they should all be real as S2 is symmetric)
		if (d(0) <= 0 || d(1) <= 0 || d(2) <= 0)
		{
			//d(0) = std::max(d(0), real(0));
			//d(1) = std::max(d(1), real(0));
			//d(2) = std::max(d(2), real(0));
			return false;
		}

		Vector3R v0 = Vector3R(V(0, 0), V(1, 0), V(2, 0));
		Vector3R v1 = Vector3R(V(0, 1), V(1, 1), V(2, 1));
		Vector3R v2 = Vector3R(V(0, 2), V(1, 2), V(2, 2));
		S = sqrt(d(0)) * Matrix3R::TensorProduct(v0, v0) +
			sqrt(d(1)) * Matrix3R::TensorProduct(v1, v1) +
			sqrt(d(2)) * Matrix3R::TensorProduct(v2, v2);
		//R = A * S.GetInverse();

		Matrix3R Sinv;
		S = (1.f / sqrt(d(0))) * Matrix3R::TensorProduct(v0, v0) +
			(1.f / sqrt(d(1))) * Matrix3R::TensorProduct(v1, v1) +
			(1.f / sqrt(d(2))) * Matrix3R::TensorProduct(v2, v2);
		R = A * Sinv;

		return true;
	}

	void PolarDecomposition::QRPolarDecomposition(const Matrix3R& F, Matrix3R& R)
	{
		// Gram-Schmidt
		Vector3R r0 = F(0);
		r0.Normalize();
		Vector3R a1 = F(1);
		Vector3R r1 = a1 - a1.Dot(r0) * r0;
		r1.Normalize();
		Vector3R r2 = r0.Cross(r1);
		r2.Normalize();
		R = Matrix3R(r0, r1, r2);
	}

}