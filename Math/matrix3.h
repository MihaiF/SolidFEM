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

#ifndef MATRIX3_H
#define MATRIX3_H

#include "Vector3.h"

using namespace Math; // not a very good idea

// TODO: finish removing float and replacing it by Real

template<typename Real>
struct Matrix3T
{

	static Matrix3T Identity() { static Matrix3T id; return id; }
	static Matrix3T Zero() { return Matrix3T(0, 0, 0, 0, 0, 0, 0, 0, 0); }

	Real m[3][3]; // TODO: linearize - column major

	static Matrix3T RotationY(float u)
	{
		float c = cos(u);
		float s = sin(u);
		Matrix3T R;
		R.m[0][0] = c;
		R.m[0][2] = s;
		R.m[2][0] = -s;
		R.m[2][2] = c;
		return R;
	}

	static Matrix3T RotationX(float u)
	{
		float c = cos(u);
		float s = sin(u);
		Matrix3T R;
		R.m[1][1] = c;
		R.m[1][2] = s;
		R.m[2][1] = -s;
		R.m[2][2] = c;
		return R;
	}

	static Matrix3T RotationZ(float u)
	{
		float c = cos(u);
		float s = sin(u);
		Matrix3T R;
		R.m[0][0] = c;
		R.m[0][1] = s;
		R.m[1][0] = -s;
		R.m[1][1] = c;
		return R;
	}

	Matrix3T GetTranspose() const
	{
		Matrix3T ret;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				ret.m[i][j] = m[j][i];
			}
		}
		return ret;
	};

	// construct matrix from column vectors
	Matrix3T(const Vector3T<Real>& v0, const Vector3T<Real>& v1, const Vector3T<Real>& v2)
	{
		for (int i = 0; i < 3; i++)
		{
			m[i][0] = v0[i];
			m[i][1] = v1[i];
			m[i][2] = v2[i];
		}
	}

	Matrix3T(Real a00, Real a01, Real a02,
		Real a10, Real a11, Real a12,
		Real a20, Real a21, Real a22)
	{
		m[0][0] = a00;
		m[0][1] = a01;
		m[0][2] = a02;
		m[1][0] = a10;
		m[1][1] = a11;
		m[1][2] = a12;
		m[2][0] = a20;
		m[2][1] = a21;
		m[2][2] = a22;
	}

	Matrix3T(const Vector3T<Real>& v)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
					m[i][j] = v[i];
				else m[i][j] = 0;
			}
		}
	}

	Matrix3T(Real d)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
					m[i][j] = d;
				else m[i][j] = 0;
			}
		}
	}

	Matrix3T() {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) {
				if (i == j) m[i][j] = 1;
				else m[i][j] = 0;
			}
	}

	Matrix3T(const Matrix3T& M) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m[i][j] = M.m[i][j];
	}

	Matrix3T& operator=(const Matrix3T& M) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m[i][j] = M.m[i][j];
		return *this;
	}

	Vector3T<Real> operator()(int i) const
	{
		Vector3T<Real> v;
		for (int j = 0; j < 3; j++)
			v[j] = m[j][i];
		return v;
	}

	Vector3T<Real> operator[](int i) const
	{
		Vector3T<Real> v;
		for (int j = 0; j < 3; j++)
			v[j] = m[i][j];
		return v;
	}

	Real operator()(int i, int j) const
	{
		return m[i][j];
	}

	Real& operator()(int i, int j)
	{
		return m[i][j];
	}

	Vector3T<Real> GetDiagonal() {
		Vector3T<Real> v;
		for (int i = 0; i < 3; i++)
			v[i] = m[i][i];
		return v;
	}

	Real Determinant()
	{
		Real a00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
		Real a01 = -m[1][0] * m[2][2] + m[2][0] * m[1][2];
		Real a02 = m[1][0] * m[2][1] - m[2][0] * m[1][1];
		return m[0][0] * a00 + m[0][1] * a01 + m[0][2] * a02;
	}

	Real Trace()
	{
		return m[0][0] + m[1][1] + m[2][2];
	}

	void Scale(float s)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m[i][j] *= s;
	}

	void Scale(double s)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m[i][j] *= (Real)s;
	}

	// returns the inverse of this matrix (without modifying it)
	Matrix3T GetInverse() const
	{
		Matrix3T inv;
		inv.m[0][0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];
		inv.m[1][0] = -m[1][0] * m[2][2] + m[2][0] * m[1][2];
		inv.m[2][0] = m[1][0] * m[2][1] - m[2][0] * m[1][1];
		inv.m[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];
		inv.m[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];
		inv.m[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];
		inv.m[0][2] = m[0][1] * m[1][2] - m[0][2] * m[1][1];
		inv.m[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];
		inv.m[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];
		Real det = m[0][0] * inv.m[0][0] + m[0][1] * inv.m[1][0] + m[0][2] * inv.m[2][0];
		inv.Scale(1.f / det);
		return inv;
	}

	static Matrix3T TensorProduct(const Vector3T<Real>& a, const Vector3T<Real>& b)
	{
		Matrix3T ret;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				ret.m[i][j] = a[i] * b[j];
			}
		}
		return ret;
	}

	static Matrix3T Skew(const Vector3T<Real>& v)
	{
		return Matrix3T(0, -v.z, v.y,
						v.z, 0, -v.x,
						-v.y, v.x, 0);
	}

	Matrix3T operator-() const
	{
		Matrix3T B;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				B.m[i][j] = -m[i][j];
		return B;
	}

	static Real DoubleContraction(const Matrix3T& A, const Matrix3T& B)
	{
		Real sum = 0;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				sum += A.m[i][j] * B.m[i][j];
			}
		}
		return sum;
	}
};

template<typename Real>
inline Vector3T<Real> vTransform(const Matrix3T<Real>& M, const Vector3T<Real>& v)
{
	Vector3T<Real> r;
	for (int i = 0; i < 3; i++) 
	{
		r[i] = 0;
		for (int j = 0; j < 3; j++)
			r[i] += M.m[i][j] * v[j];
	}
	return r;
}

template<typename Real>
inline Vector3T<Real> operator*(const Matrix3T<Real>& M, const Vector3T<Real>& v)
{
	return vTransform(M, v);
}

template<typename Real>
inline Matrix3T<Real> operator*(float s, const Matrix3T<Real>& M)
{
	Matrix3T<Real> ret = M;
	ret.Scale(s);
	return ret;
}

template<typename Real>
inline Matrix3T<Real> operator*(double s, const Matrix3T<Real>& M)
{
	Matrix3T<Real> ret = M;
	ret.Scale(s);
	return ret;
}

template<typename Real>
inline Matrix3T<Real> operator*(const Matrix3T<Real>& A, const Matrix3T<Real>& B)
{
	Matrix3T<Real> C;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			C.m[i][j] = 0;
			for (int k = 0; k < 3; k++)
				C.m[i][j] += A.m[i][k] * B.m[k][j];
		}
	}
	return C;
}

template<typename Real>
inline Matrix3T<Real> mTranspose(const Matrix3T<Real>& A)
{
	Matrix3T<Real> B;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			B.m[i][j] = A.m[j][i];
	return B;
}

template<typename Real>
inline Matrix3T<Real> operator!(const Matrix3T<Real>& A)
{
	return mTranspose(A);
}

template<typename Real>
inline Matrix3T<Real> operator+(const Matrix3T<Real>& A, const Matrix3T<Real>& B) 
{
	Matrix3T<Real> C;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			C.m[i][j] = A.m[i][j] + B.m[i][j];
	return C;
}

template<typename Real>
inline Matrix3T<Real> operator-(const Matrix3T<Real>& A, const Matrix3T<Real>& B)
{
	Matrix3T<Real> C;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			C.m[i][j] = A.m[i][j] - B.m[i][j];
	return C;
}

typedef Matrix3T<float> Matrix3;

// TODO: convert to templated version (or drop)
Matrix3 mDiagonal(float a, float b, float c);
Vector3 vTransformTr(Matrix3 M, Vector3 v);
Vector3 operator*(Vector3, Matrix3 M);
Matrix3 operator~(Vector3 v); // skew symmetric tensor

void GetTransform(float T[16], Matrix3 R, Vector3 X);

inline Matrix3 mDiagonal(float a, float b, float c)
{
	Matrix3 M;
	M.m[0][0]=a;
	M.m[1][1]=b;
	M.m[2][2]=c;
	//presupun ca celelate elem sunt 0 (deocamdata)
	return M;
}

inline Vector3 vTransformTr(Matrix3 M, Vector3 v){
	Vector3 r;
	for(int i=0;i<3;i++){
		r[i]=0;
		for(int j=0;j<3;j++)
			r[i]+=M.m[j][i]*v[j];
	}
	return r;
}

inline Vector3 operator*(Vector3 v, Matrix3 M){
	return vTransformTr(M,v);
}

inline Matrix3 operator~(Vector3 v){
	Matrix3 M;
	M.m[0][1] = -v.Z();
	M.m[1][0] = v.Z();
	M.m[0][2] = v.Y();
	M.m[2][0] = -v.Y();
	M.m[1][2] = -v.X();
	M.m[2][1] = v.X();
	return M;
}

inline void GetTransform(float T[16], Matrix3 R, Vector3 X)
{
  T[0] = R.m[0][0];
  T[1] = R.m[1][0];
  T[2] = R.m[2][0];
  T[3] = 0;
  T[4] = R.m[0][1];
  T[5] = R.m[1][1];
  T[6] = R.m[2][1];
  T[7] = 0;
  T[8] = R.m[0][2];
  T[9] = R.m[1][2];
  T[10] = R.m[2][2];
  T[11] = 0;
  T[12] = X.X();
  T[13] = X.Y();
  T[14] = X.Z();
  T[15] = 1;	
}

#endif // MATRIX3_H