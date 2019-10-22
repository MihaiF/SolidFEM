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

#ifndef FEM_DATA_STRUCTURE_H
#define FEM_DATA_STRUCTURE_H	
#include <Engine/Types.h>
#include <Math/Vector3.h>
#include <Math/Matrix3.h>
#include <vector>

#if !defined(_DEBUG) && defined(USE_MKL)
	#define EIGEN_USE_MKL_ALL
	#include <Eigen/PardisoSupport>
#endif
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#define USE_DOUBLE_FOR_FEM

namespace FEM_SYSTEM
{
#ifdef USE_DOUBLE_FOR_FEM
	typedef double real;
#else
	typedef float real;
#endif

	typedef Matrix3T<real> Matrix3R;
	typedef Vector3T<real> Vector3R;

	typedef Eigen::SparseMatrix<real> SparseMatrix;

	typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
	typedef Eigen::Matrix<real, Eigen::Dynamic, 1> EigenVector;

	enum SimulationType
	{
		ST_STATIC,
		ST_QUASI_STATIC,
		ST_EXPLICIT,
		ST_IMPLICIT,
	};

	enum MethodType
	{
		MT_LINEAR_ELASTICITY,
		MT_COROTATIONAL_ELASTICITY,
		MT_INCOMPRESSIBLE_LINEAR_ELASTICITY,
		MT_INCOMPRESSIBLE_COROTATIONAL_ELASTICITY,
		MT_NONLINEAR_ELASTICITY,
		MT_CONSTRAINT_LINEAR,
		MT_ANALYTIC_CANTILEVER_NONLINEAR_ELASTICITY,
		MT_ANALYTIC_CANTILEVER_LINEAR_ELASTICITY,
		MT_VEGA
	};

	struct FemConfig
	{
		real mYoungsModulus = 66000;
		real mPoissonRatio = 0.3f;
		real mDampingYoungsModulus = 0;
		real mDampingPoissonRatio = 0;
		MethodType mType = MT_LINEAR_ELASTICITY;
		int mOrder = 1;
		real mDensity = 1070.f;
		real mGravity = -9.8f;
		SimulationType mSimType = ST_STATIC;
		int mNumSubsteps = 1;
		real mForceApplicationStep = 0.1f;
		void* mCustomConfig = nullptr;
	};

	struct Node
	{
		Node() : invMass(1.f) {}

		Vector3R pos, vel, force;
		Vector3R pos0; // undeformed position
		real invMass;
	};

	struct Tetrahedron
	{
		uint16 i[4]; // global node indices
		Matrix3R X, Xtr; // undeformed shape matrix inverse (and its transpose)
		Matrix3R R; // used by SubStepLinearImplicitCG
		real vol; // volume
		// linear FEM
		Matrix3R K[4][4]; // stiffness matrix
		Matrix3R Hn[4], Hs[4]; // fixed Jacobian terms
		// FVM
		Vector3R NA[4]; // undeformed oriented face areas
		//Vector3R b[4];
		//Matrix3R Bm;

		Tetrahedron() { }
		Tetrahedron(uint16 i0, uint16 i1, uint16 i2, uint16 i3)
		{
			i[0] = i0;
			i[1] = i1;
			i[2] = i2;
			i[3] = i3;
		}
	};

	typedef Eigen::Matrix<real, 6, 1> Vector6;

	// A symmetric tensor (e.g. strain or stress) composed of 3 normal components and 3 shear components
	struct SymTensor
	{
		Vector3R diag, shear;
	};

	// Helper structure used by the velocity solver to pre-cache data
	struct TetInfo
	{
		Eigen::Matrix<real, 6, 12> J;
		Eigen::Matrix<real, 6, 6> S;
		Eigen::Matrix<real, 6, 1> c;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	inline SymTensor Tensor2Vector(const Matrix3R& M)
	{
		SymTensor ret;
		ret.diag.Set(M.m[0][0], M.m[1][1], M.m[2][2]);
		ret.shear.Set(M.m[0][1], M.m[1][2], M.m[2][0]);
		return ret;
	};

	inline Matrix3R Vector2Tensor(const SymTensor& st)
	{
		return Matrix3R(st.diag.X(), st.shear.X(), st.shear.Z(),
			st.shear.X(), st.diag.Y(), st.shear.Y(),
			st.shear.Z(), st.shear.Y(), st.diag.Z());
	};

	template<typename MATRIX>
	inline MATRIX Matrix3ToEigen(const Matrix3R& R)
	{
		EigenMatrix Reig(3, 3);
		Reig(0, 0) = R(0, 0);
		Reig(0, 1) = R(0, 1);
		Reig(0, 2) = R(0, 2);
		Reig(1, 0) = R(1, 0);
		Reig(1, 1) = R(1, 1);
		Reig(1, 2) = R(1, 2);
		Reig(2, 0) = R(2, 0);
		Reig(2, 1) = R(2, 1);
		Reig(2, 2) = R(2, 2);
		return Reig;
	}

	//TO BE MOVED IN A CLASS!!!

	struct MeshInterp
	{
		int id1, id2, id3;
		float u, v, w;//barycentric
		Vector3 point;
		Vector3 original;
		float distance;
		int tet_index;
		int vert_index;
	};


inline Vector3 ClosestPtPointTriangle(const Vector3& p,
								   const Vector3& a,
								   const Vector3& b,
								   const Vector3& c)
{
	// Check if P in vertex region outside A
	Vector3 ab = b - a;
	Vector3 ac = c - a;
	Vector3 ap = p - a;
	float d1 = ab.Dot(ap);
	float d2 = ac.Dot(ap);
	if (d1 <= 0.0f && d2 <= 0.0f) return a; // barycentric coordinates (1,0,0)
	// Check if P in vertex region outside B
	Vector3 bp = p - b;
	float d3 = ab.Dot(bp);
	float d4 = ac.Dot(bp);
	if (d3 >= 0.0f && d4 <= d3) return b; // barycentric coordinates (0,1,0)
	// Check if P in edge region of AB, if so return projection of P onto AB
	float vc = d1*d4 - d3*d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
		float v = d1 / (d1 - d3);
		return a + v * ab; // barycentric coordinates (1-v,v,0)
	}
	// Check if P in vertex region outside C
	Vector3 cp = p - c;
	float d5 = ab.Dot(cp);
	float d6 = ac.Dot(cp);
	if (d6 >= 0.0f && d5 <= d6) return c; // barycentric coordinates (0,0,1)
	// Check if P in edge region of AC, if so return projection of P onto AC
	float vb = d5*d2 - d1*d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
		float w = d2 / (d2 - d6);
		return a + w * ac; // barycentric coordinates (1-w,0,w)
	}
	// Check if P in edge region of BC, if so return projection of P onto BC
	float va = d3*d6 - d5*d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
		float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return b + w * (c - b); // barycentric coordinates (0,1-w,w)
	}
	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	float denom = 1.0f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;
	return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
}

inline int PointOutsideOfPlane(const Vector3& p,
							const Vector3& a,
							const Vector3& b,
							const Vector3& c,
							const Vector3& d)
{
	float signp = (p - a).Dot((b - a).Cross(c - a)); // [AP AB AC]
	float signd = (d - a).Dot((b - a).Cross(c - a)); // [AD AB AC]
	// Points on opposite sides if expression signs are opposite
	return signp * signd < 0.0f;
}

	inline int PointOutsideOfPlane(const Vector3& p, 
							const Vector3& a, 
							const Vector3& b,
							const Vector3& c)
	{
		return (p - a).Dot((b - a).Cross(c - a)) >= 0.0f; // [AP AB AC] >= 0
	}

	inline float ClosestPtPointTetrahedron(const Vector3& p, 
									const Vector3& a,
									const Vector3& b,
									const Vector3& c,
									const Vector3& d,
									int& id1,
									int& id2,
									int& id3,
									Vector3& point)
{
	// Start out assuming point inside all halfspaces, so closest to itself
	int count = 0;
	float bestSqDist = FLT_MAX;
	// If point outside face abc then compute closest point on abc
	if (PointOutsideOfPlane(p, a, b, c))
	{
		Vector3 q = ClosestPtPointTriangle(p, a, b, c);
		float sqDist = (q - p).Dot(q - p);
		// Update best closest point if (squared) distance is less than current best
		if (sqDist < bestSqDist){ bestSqDist = sqDist; point = q; id1 = 0; id2 = 1; id3 = 2; }
	}
	else
	{
		count++;
	}
	// Repeat test for face acd
	if (PointOutsideOfPlane(p, a, c, d)) {
		Vector3 q = ClosestPtPointTriangle(p, a, c, d);
		float sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist){ bestSqDist = sqDist; point = q; id1 = 0; id2 = 2; id3 = 3; }
	}
	else
	{
		count++;
	}
	// Repeat test for face adb
	if (PointOutsideOfPlane(p, a, d, b)) {
		Vector3 q = ClosestPtPointTriangle(p, a, d, b);
		float sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist){ bestSqDist = sqDist; point = q; id1 = 0; id2 = 3; id3 = 1; }
	}
	else
	{
		count++;
	}
	// Repeat test for face bdc
	if (PointOutsideOfPlane(p, b, d, c)) {
		Vector3 q = ClosestPtPointTriangle(p, b, d, c);
		float sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist) { bestSqDist = sqDist; point = q;  id1 = 1; id2 = 3; id3 = 2; }
	}
	else
	{
		count++;
	}

	//point inside tetrahedron
	if (count == 4)
	{
		Vector3 q = ClosestPtPointTriangle(p, a, b, c);
		float sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist){ bestSqDist = sqDist; point = q; id1 = 0; id2 = 1; id3 = 2; }

		q = ClosestPtPointTriangle(p, a, c, d);
		sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist){ bestSqDist = sqDist; point = q; id1 = 0; id2 = 2; id3 = 3; }

		q = ClosestPtPointTriangle(p, a, d, b);
		sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist){ bestSqDist = sqDist; point = q; id1 = 0; id2 = 3; id3 = 1; }

		q = ClosestPtPointTriangle(p, b, d, c);
		sqDist = (q - p).Dot(q - p);
		if (sqDist < bestSqDist) { bestSqDist = sqDist; point = q;  id1 = 1; id2 = 3; id3 = 2; }

	}

	return bestSqDist;
 }


 inline void BarycentricCoordinates(const Vector3& p,
							const Vector3& a,
							const Vector3& b,
							const Vector3& c,
							float& u,
							float& v,
							float& w)
{
		Vector3 v0 = c - a;
		Vector3 v1 = b - a;
		Vector3 v2 = p - a;

		// Compute dot products
		float dot00 = v0.Dot(v0);
		float dot01 = v0.Dot(v1);
		float dot02 = v0.Dot(v2);
		float dot11 = v1.Dot(v1);
		float dot12 = v1.Dot(v2);

		// Compute barycentric coordinates
		float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
		u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		v = (dot00 * dot12 - dot01 * dot02) * invDenom;
		w = 1.0f - u - v;
}

} // namespace FEM_SYSTEM

#endif