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
#include <iostream>


#if !defined(_DEBUG) && defined(USE_MKL)
	#define EIGEN_USE_MKL_ALL
	#include <Eigen/PardisoSupport>
#endif
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#define USE_DOUBLE_FOR_FEM

// put them in the global namespace until we create a new namespace
#ifdef USE_DOUBLE_FOR_FEM
typedef double real;
#else
typedef float real;
#endif

typedef Eigen::SparseMatrix<real> SparseMatrix;

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
typedef Eigen::Matrix<real, Eigen::Dynamic, 1> EigenVector;

typedef Eigen::Matrix<real, 3, 3> EigenMatrix3;

namespace FEM_SYSTEM
{
	typedef Matrix3T<real> Matrix3R;
	typedef Vector3T<real> Vector3R;

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
		MT_INCOMPRESSIBLE_NONLINEAR_ELASTICITY,
		MT_CONSTRAINT_LINEAR,
		MT_ANALYTIC_CANTILEVER_NONLINEAR_ELASTICITY,
		MT_ANALYTIC_CANTILEVER_LINEAR_ELASTICITY,
	};

	enum MaterialModelType
	{
		MMT_LINEAR,
		MMT_STVK,
		MMT_COROTATIONAL,
		MMT_NEO_HOOKEAN,
		MMT_NEO_HOOKEAN_OGDEN,
		MMT_DISTORTIONAL_LINEAR,
		MMT_DISTORTIONAL_MOONEY, 
		MMT_DISTORTIONAL_NH2,
		MMT_DISTORTIONAL_NH3,
		MMT_DISTORTIONAL_BW,
		MMT_DISTORTIONAL_OGDEN,
		MMT_DISTORTIONAL_NH6,
		MMT_DISTORTIONAL_NH7,
	};

	enum NonlinearSolverType
	{
		NST_NEWTON,
		NST_NEWTON_LS,
		NST_NEWTON_CG,
		NST_NONLINEAR_CG,
		NST_GRADIENT_DESCENT,
		NST_STEEPEST_DESCENT,
		NST_NLOPT,
		NST_GSL,
		NST_MCL,
		NST_MIXED_NEWTON,
		NST_MIXED_DUAL_ASCENT,
	};

	enum // some integer constants
	{
		NUM_STRESS_COMPONENTS = 6,
		NUM_POS_COMPONENTS = 3,
		NUM_BARYCENTRIC_COMPONENTS = 4,
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
		real mForceApplicationStep = 0.1;
		void* mCustomConfig = nullptr;
		MaterialModelType mMaterial = MMT_LINEAR;
		uint32 mOuterIterations = 10;
		uint32 mInnerIterations = 100;
		bool mHasCollisions = false;
		real mInvBulkModulus = -1;
		real mContactStiffness = 20000;
		uint32 mNumFixed = 0;
		real mAbsNewtonRsidualThreshold = 0.1;
		bool mVerbose = true;
		real mAppliedPressure = 0; // for loading purposes only for now
		bool mZUpAxis = false;
	};

	struct Node
	{
		Node() : invMass(1.f) {}

		Vector3R pos, vel, force;
		Vector3R pos0; // undeformed position
		real invMass;
	};

	struct Tet
	{
		uint32 idx[4];
	};

	struct Tetrahedron
	{
		uint32 i[4]; // global node indices
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
		Tetrahedron(uint32 i0, uint32 i1, uint32 i2, uint32 i3)
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

	typedef std::vector<Vector3R> Vector3Array;

	inline Eigen::Map<EigenVector> GetEigenVector(Vector3Array& arr, uint32 start = 0)
	{
		return Eigen::Map<EigenVector>((real*)&arr[start], (arr.size() - start) * 3, 1);
	}

	inline Vector3Array GetStdVector(const EigenVector& vec)
	{
		size_t size = vec.size() / 3;
		Vector3R* ptr = (Vector3R*)vec.data();
		return Vector3Array(ptr, ptr + size);
	}

} // namespace FEM_SYSTEM

#endif