#ifndef FEM_PHYSICS_LINEAR_ELASTICITY_H
#define FEM_PHYSICS_LINEAR_ELASTICITY_H

#define MEASURE_TIME_P(...)
#define PROFILE_SCOPE(...)

#include "FemPhysicsBase.h"
//#include "TetrahedralMesh.h"
#include <memory>

namespace FEM_SYSTEM
{
	typedef std::vector<Vector3R> Vector3Array;

	inline Eigen::Map<EigenVector> GetEigenVector(Vector3Array& arr)
	{
		return Eigen::Map<EigenVector>((real*)&arr[0], arr.size() * 3, 1);
	}

	inline Vector3Array GetStdVector(EigenVector& vec)
	{
		size_t size = vec.size() / 3;
		Vector3R* ptr = (Vector3R*)vec.data();
		return Vector3Array(ptr, ptr + size);
	}

	class FemPhysicsLinear : public FemPhysicsBase
	{
	public:
		enum // some integer constants
		{
			NUM_STRESS_COMPONENTS = 6,
			NUM_POS_COMPONENTS = 3,
			NUM_BARYCENTRIC_COMPONENTS = 4,
		};

		struct BarycentricJacobian
		{
			Vector3R y[4]; // TODO: store only 3 as y0 = -y1 - y2 - y3
		};

	public:
		FemPhysicsLinear(std::vector<Tetrahedron>& tetrahedra,
			std::vector<Node>& allNodes, const FemConfig& config);

		uint32 GetNumNodes() const override { return mTetMesh->GetNumNodes(); }
		uint32 GetNumFreeNodes() const { ASSERT(mNumBCs < GetNumNodes()); return GetNumNodes() - mNumBCs; }
		uint32 GetNumLocalNodes() const { return mTetMesh->GetNumNodesPerElement(); }
		uint32 GetNumElements() const { return  mTetMesh->GetNumElements(); }
		bool IsNodeFixed(uint32 i) const override { ASSERT(i < GetNumNodes()); return mReshuffleMap[i] < mNumBCs; }
		Vector3R GetDeformedPosition(uint32 i) const override;
		void SetBoundaryConditionsSurface(const std::vector<uint32>& triangleList, const std::vector<uint32>& elemList, real pressure) override;
		
		TetrahedralMesh<uint32>* GetTetMesh() { return mTetMesh.get(); }

		void SetExternalForces(const std::vector<Vector3R>& forces) override;
		real GetTotalVolume() const;

	protected:		
		void ComputeLocalMassMatrixBB2(real density, uint32 numLocalNodes, EigenMatrix& Mlocal);
		void ComputeLocalMassMatrixBB(real density, uint32 numLocalNodes, EigenMatrix& Mlocal);
		void AssembleMassMatrix();
		void CreateMeshAndDofs(std::vector<Tetrahedron>& tetrahedra, std::vector<Node>& nodes);

		void ComputeBarycentricJacobian(uint32 i, Vector3R y[4]);
		void ComputeStrainJacobian(uint32 i, Matrix3R Bn[4], Matrix3R Bs[4]);
		void ComputeBodyForces(Vector3Array& f);

		void ComputeDeformationGradient(uint32 e, Matrix3R& F);

		void ComputeLocalStiffnessMatrixBB(uint32 elementidx, real mu, real lambda, EigenMatrix& Klocal);

		void ComputeTractionForces();
		EigenVector ComputeLoadingForces();

		bool CheckForInversion();

	protected:
		uint32 mOrder;
		std::unique_ptr<TetrahedralMesh<uint32>> mTetMesh; // the topological data

		std::vector<real> mElementVolumes; // signed volumes of elements
		std::vector<Vector3R> mDeformedPositions; // nodal spatial positions
		std::vector<Vector3R> mReferencePositions; // nodal material positions
		std::vector<Vector3R> mVelocities; // nodal spatial velocities
		std::vector<BarycentricJacobian> mBarycentricJacobians; // the Jacobians of the sub-parametric mappings

		uint32 mNumBCs;
		SparseMatrix mMassMatrix; // mass matrix

		std::vector<uint32> mReshuffleMapInv; // mapping from shuffled nodes to original ones
		std::vector<uint32> mReshuffleMap; // mapping from original nodes to shuffled ones

		Vector3Array mBodyForces;
		Vector3Array mExternalForces;

		// pressure BCs
		real mAppliedPressure;
		bool mUseImplicitPressureForces = false; // compute a stiffness matrix due to the variation of the normal
		Vector3Array mTractionForces; // resulting nodal forces from Neuman BCs
		std::vector<uint32> mTractionSurface; // triangle list of Neumann boundary
		EigenMatrix mTractionStiffnessMatrix;

	private:
		bool mUseLumpedMass;

		friend class FemTester;
	};

	inline Vector3R FemPhysicsLinear::GetDeformedPosition(uint32 i) const
	{
		uint32 idx = mReshuffleMap[i];
		return idx >= mNumBCs ? mDeformedPositions[idx - mNumBCs] : mReferencePositions[idx];
	}

	// the combinations formula
	inline uint32 Combinations(int n, int k)
	{
		if (n < 0 || k < 0)
			return 0;
		// TODO: use a look up table
		uint32 denom = 1;
		uint32 nom = 1;
		for (int i = 1; i <= k; i++)
		{
			denom *= i;
			nom *= (n - i + 1);
		}
		return nom / denom;
	}

	struct MultiIndex
	{
		int mIndices[4];

		MultiIndex(int* ptr) { memcpy(mIndices, ptr, 4 * sizeof(int)); }
		int operator[](int idx) const { ASSERT(idx < 4); return mIndices[idx]; }
		int& operator[](int idx) { ASSERT(idx < 4); return mIndices[idx]; }
		int GetSum() const { return mIndices[0] + mIndices[1] + mIndices[2] + mIndices[3]; }
		bool IsNegative() const { return mIndices[0] < 0 || mIndices[1] < 0 || mIndices[2] < 0 || mIndices[3] < 0; }
		bool Decrement(uint32 idx)
		{
			ASSERT(idx < 4);
			mIndices[idx]--;
			return mIndices[idx] >= 0;
		}
	};

	// the G function from [Weber]
	inline real ComputeMultiIndexSumFactor(uint32 order, MultiIndex multiIndexI, MultiIndex multiIndexJ)
	{
		uint32 c1 = Combinations(multiIndexI[0] + multiIndexJ[0], multiIndexI[0]);
		uint32 c2 = Combinations(multiIndexI[1] + multiIndexJ[1], multiIndexI[1]);
		uint32 c3 = Combinations(multiIndexI[2] + multiIndexJ[2], multiIndexI[2]);
		uint32 c4 = Combinations(multiIndexI[3] + multiIndexJ[3], multiIndexI[3]);
		real denom = (real)Combinations(2 * order, order);
		return c1 * c2 * c3 * c4 / denom;
	}

	inline int Binom(int n, int k)
	{
#ifndef _UNIT_TESTING
		ASSERT(n >= 0);
		ASSERT(k >= 0);
		ASSERT(n >= k);
#endif // !_UNIT_TESTING
		if (n < 0)
		{
			throw std::invalid_argument("Binom is undefined for n < 0");
		}
		if (k < 0)
		{
			throw std::invalid_argument("Binom is undefined for k < 0");
		}
		if (n < k)
		{
			throw std::invalid_argument("Binom is undefined for n < k");
		}

		int val = 1;

		// Utilize that C(n, k) = C(n, n-k) 
		if (k > n - k)
			k = n - k;

		for (int i = 0; i < k; ++i)
		{
			val *= (n - i);
			val /= (i + 1);
		}

		return val;
	}

	inline real G(int ijkl[4], int mnop[4], int order)
	{
		ASSERT(ijkl[0] + ijkl[1] + ijkl[2] + ijkl[3] == order);
		ASSERT(mnop[0] + mnop[1] + mnop[2] + mnop[3] == order);

		if (ijkl[0] < 0 || ijkl[1] < 0 || ijkl[2] < 0 || ijkl[3] < 0 ||
			mnop[0] < 0 || mnop[1] < 0 || mnop[2] < 0 || mnop[3] < 0)
			return 0;

		int b1 = Binom(ijkl[0] + mnop[0], ijkl[0]);
		int b2 = Binom(ijkl[1] + mnop[1], ijkl[1]);
		int b3 = Binom(ijkl[2] + mnop[2], ijkl[2]);
		int b4 = Binom(ijkl[3] + mnop[3], ijkl[3]);
		real b1234 = real(b1 * b2 * b3 * b4);
		int bOrder = Binom(2 * order, order);
		return real(b1234 / bOrder);
	}

	inline void PrintMatrix(const EigenMatrix& M)
	{
		for (int i = 0; i < M.rows(); i++)
		{
			for (int j = 0; j < M.cols(); j++)
			{
				Printf("%g\t\t", M(i, j));
			}
			Printf("\n");
		}
	}

	inline Matrix3R Symm(const Vector3R& v)
	{
		return Matrix3R(0, v.Z(), v.Y(),
			v.Z(), 0, v.X(),
			v.Y(), v.X(), 0);
	}

	inline void FemPhysicsLinear::SetExternalForces(const std::vector<Vector3R>& forces)
	{
		mExternalForces.resize(GetNumFreeNodes());
		for (uint32 i = 0; i < forces.size(); i++)
		{
			uint32 idx = mReshuffleMap[i];
			if (idx >= mNumBCs)
				mExternalForces[idx - mNumBCs] = forces[i];
		}
	}

} // namespace FEM_SYSTEM

#endif // FEM_PHYSICS_LINEAR_ELASTICITY_H
