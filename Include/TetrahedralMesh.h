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

#ifndef TETRAHEDRAL_MESH_H
#define TETRAHEDRAL_MESH_H

// C++ references
#include <tuple>

// Custom references
#include "StridedVector.h"
#include "ITetrahedralMesh.h"
#include <Engine/Types.h>
#include <Engine/Utils.h>

namespace FEM_SYSTEM
{
	template<class IndexType>
	class TetrahedralMesh : public ITetrahedralMesh
	{
	public:

		// Typedefs
		typedef std::tuple<IndexType, int, int> InterpolationTuple;
		typedef std::tuple<IndexType, IndexType, IndexType> tuple3IndexType;
		typedef std::tuple<IndexType, IndexType, IndexType, IndexType> tuple4IndexType;
		typedef IndexType mIndexType;

		// 'uint8' ensures that we support an order up to 255
		// -> I,J,K,L are always \in [0, order]
		//typedef uint8 IJKLType;
		typedef int IJKLType; // TODO need to fix MultiIndex before anything but int can be used

		struct CompressedElement
		{
			IndexType c1, c2, c3, c4;
			IndexType e1;
			IndexType e2;
			IndexType e3;
			IndexType e4;
			IndexType e5;
			IndexType e6;

			IndexType f1;
			IndexType f2;
			IndexType f3;
			IndexType f4;

			IndexType b1;

			// TODO: use std::bitset instead of raw bits, ex std::bitset<18>? Depends, are there any performance advantages?
			// -> Need performance tests to evaluate that.
			uint8 edge_reverse_mask : 6;
			uint8 face_bucketflip_mask : 4;
			uint8 face_rotation_mask : 8;

			CompressedElement()
			{
				// TODO should we initialize all of the gidx values? with 0?
				edge_reverse_mask = 0;
				face_rotation_mask = 0;
				face_bucketflip_mask = 0;
			}

			// clidx \in [0,3]
			IndexType GetCornerNode(int clidx)
			{
				ASSERT(clidx >= 0);
				ASSERT(clidx < 4);

				IndexType *node = &this->c1;
				return *(IndexType*)((char*)node + sizeof(IndexType)*clidx);
			}

			// clidx \in [0,3]
			void SetCornerNode(int clidx, IndexType gidx)
			{
				ASSERT(clidx >= 0);
				ASSERT(clidx < 4);

				IndexType *node = &this->c1;
				*(IndexType*)((char*)node + sizeof(IndexType)*clidx) = gidx;
			}

			// elidx \in [0,5]
			void GetEdgeNode(int elidx, IndexType* gidx, bool* reverse)
			{
				ASSERT(elidx >= 0);
				ASSERT(elidx < 6);

				IndexType *p_gidx = &this->e1;
				p_gidx = (IndexType*)((char*)p_gidx + sizeof(IndexType)*elidx);
				*gidx = *p_gidx;

				*reverse = (edge_reverse_mask >> elidx) & (uint8)1;
			}

			// elidx \in [0,5]
			void SetEdgeNode(int elidx, IndexType gidx, bool reverse)
			{
				ASSERT(elidx >= 0);
				ASSERT(elidx < 6);

				IndexType *p_gidx = &this->e1;
				p_gidx = (IndexType*)((char*)p_gidx + sizeof(IndexType)*elidx);
				*p_gidx = gidx;

				reverse = !!reverse; // Ensure that reverse \in [0,1], really needed?
				edge_reverse_mask = edge_reverse_mask & ~(1 << elidx) | ((uint8)reverse << elidx);
			}

			// flidx \in [0,3]
			void GetFaceNode(int flidx, IndexType* gidx, int* rotations, bool* bucketflip)
			{
				ASSERT(flidx >= 0);
				ASSERT(flidx < 4);

				IndexType *p_gidx = &this->f1;
				p_gidx = (IndexType*)((char*)p_gidx + sizeof(IndexType)*flidx);
				*gidx = *p_gidx;

				*rotations = (face_rotation_mask >> (flidx * 2)) & (uint8)3;
				*bucketflip = (face_bucketflip_mask >> flidx) & (uint8)1;
			}

			// flidx \in [0,3]
			// rotations \in [0,2]
			void SetFaceNode(int flidx, IndexType gidx, int rotations, bool bucketflip)
			{
				ASSERT(flidx >= 0);
				ASSERT(flidx < 4);

				ASSERT(rotations >= 0);
				ASSERT(rotations < 3);

				IndexType *p_gidx = &this->f1;
				p_gidx = (IndexType*)((char*)p_gidx + sizeof(IndexType)*flidx);
				*p_gidx = gidx;

				face_rotation_mask = (face_rotation_mask & ~(3 << (flidx * 2))) | ((uint8)rotations << (flidx * 2));
				bucketflip = !!bucketflip; // Ensure that bucketflip \in [0,1], really needed?
				face_bucketflip_mask = face_bucketflip_mask & ~(1 << flidx) | ((uint8)bucketflip << flidx);
			}

			IndexType GetBodyNode()
			{
				return b1;
			}

			void SetBodyNode(IndexType gidx)
			{
				b1 = gidx;
			}
		};

		TetrahedralMesh();
		~TetrahedralMesh();

		void InstanceWith(int order, int noPoints, StridedVector<uint32> connectivity, int noElements);

		int GetGlobalIndex(int eidx, int nidx);

		int GetNumNodes() const;
		int GetNumElements() const;
		int GetNumNodesPerElement() const;
		int GetNumLinearNodes();

		int GetOrder() const;
		IJKLType* GetIJKL(int lidx);

		// <gidx, eidx, lidx>
		std::vector<InterpolationTuple> mInterpolationscheme;
		IJKLType* mIJKL;

		// NOTE part of experimental feature of HO_RENDERING
		int GetGlobalIndicesForSurfaceTriangle(const IndexType c1, const IndexType c2, const IndexType c3, IndexType *out, int* IJKs);

		// TODO copy pasted from naive implementation
		// Returns 0 if order < 0
		// Note that this is equivalent to Binom(q+3,3), without hasling with factorials
		static int GetNodesPerElement(int order)
		{
			int v = 0;
			for (int i = 0; i < order + 1; i++)
			{
				v += (i + 1) * (i + 2);
			}
			v = v / 2;
			return v;
		}

		int GetNodesPerEdge(int order) const
		{
			return std::max(order - 1, 0);
		}

		int GetNodesPerFace(int order) const
		{
			if (order > 0)
			{
				return std::max(((order - 2)*(order - 1)) / 2, 0);
			}
			return 0;
		}

		static int GetNodesPerBody(int order)
		{
			return GetNodesPerElement(order - 4);
		}
		// end of copy paste

		bool GetIsInstanced() { return mIsInstanced; }

	private:

		std::vector<CompressedElement> mConnectivityList;
		int mOrgNoPoints;
		int mNodesInMesh;
		int mNodesPrElement;
		int mElementsInMesh;
		int mOrder;
		bool mIsInstanced = false;

		IndexType GetGlobalIndexCorner(int eidx, int nidx);
		IndexType GetGlobalIndexEdge(int eidx, int nidx);
		IndexType GetGlobalIndexFace(int eidx, int nidx);
		IndexType GetGlobalIndexBody(int eidx, int nidx);

		void BuildHigherOrderMesh(int noPoints, StridedVector<uint32> connectivity);
		int GetRotatedFaceIdx(int idx, int order, int rotations);
		int GetBucketflippedFaceIdx(int idx, int order);
		// Including 'from' and 'to'
		int Sum(int from, int to);
	};
}



// C++ references
#include <vector>
#include <algorithm>
#include <string>

// Custom references
#include "Engine/Utils.h"

// Macros
#define MATRIX_INDEX(i, j, width) i * width + j

namespace FEM_SYSTEM
{

	template<class IndexType>
	TetrahedralMesh<IndexType>::TetrahedralMesh()
	{

	}

	template<class IndexType>
	TetrahedralMesh<IndexType>::~TetrahedralMesh()
	{
		if (mIsInstanced)
		{
			delete[] mIJKL;
		}
	}

	template<class IndexType>
	void TetrahedralMesh<IndexType>::InstanceWith(int order, int noPoints, StridedVector<uint32> connectivity, int noElements)
	{
		// TODO assert that the selected IJKLType is big enough to handle the given order
		// TODO either assert that this mesh has already been instanced, or cleanup before each instance
		ASSERT(order > 0);
		ASSERT(noPoints >= 4);
		ASSERT(noElements > 0);
		ASSERT(mIsInstanced == false);
		// TODO assert that connectivity is not null, how to?

		this->mOrder = order;
		this->mElementsInMesh = noElements;
		this->mNodesPrElement = GetNodesPerElement(order);

		this->mConnectivityList.resize(noElements);
		this->mIJKL = new IJKLType[mNodesPrElement * 4]{ 0 };
		this->mIsInstanced = true;

		BuildHigherOrderMesh(noPoints, connectivity);
		this->mNodesInMesh = mInterpolationscheme.size();
		this->mOrgNoPoints = noPoints;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetGlobalIndex(int eidx, int nidx)
	{
		ASSERT(eidx >= 0);
		ASSERT(eidx < mElementsInMesh);
		ASSERT(nidx >= 0);
		ASSERT(nidx < mNodesPrElement);

		if (nidx >= 0 && eidx >= 0 && eidx < mElementsInMesh)
		{
			int pre_width = 4;
			if (nidx < pre_width && mOrder > 0)
			{
				return GetGlobalIndexCorner(eidx, nidx);
			}
			pre_width += GetNodesPerEdge(mOrder) * 6;
			if (nidx < pre_width && mOrder > 1)
			{
				return GetGlobalIndexEdge(eidx, nidx);
			}
			pre_width += GetNodesPerFace(mOrder) * 4;
			if (nidx < pre_width && mOrder > 2)
			{
				return GetGlobalIndexFace(eidx, nidx);
			}
			pre_width += GetNodesPerBody(mOrder);
			if (nidx < pre_width && mOrder > 3)
			{
				return GetGlobalIndexBody(eidx, nidx);
			}
		}

		// The given element or local index does not exist in this mesh.
		return -1;
	}

	template<class IndexType>
	IndexType TetrahedralMesh<IndexType>::GetGlobalIndexCorner(int eidx, int nidx)
	{
		return mConnectivityList[eidx].GetCornerNode(nidx);
	}

	template<class IndexType>
	IndexType TetrahedralMesh<IndexType>::GetGlobalIndexEdge(int eidx, int nidx)
	{
		std::div_t split = std::div(nidx - 4, GetNodesPerEdge(mOrder));
		// TODO wont work if the edge_code is 0, ie if an edge node has the global idx of 0, but
		// we should be able to assume that this never occurs, right?

		IndexType e_gidx; bool edge_reverse;
		mConnectivityList[eidx].GetEdgeNode(split.quot, &e_gidx, &edge_reverse);

		if (edge_reverse)
		{
			split.rem *= -1;
		}
		return e_gidx + split.rem;
	}

	template<class IndexType>
	IndexType TetrahedralMesh<IndexType>::GetGlobalIndexFace(int eidx, int nidx)
	{
		std::div_t split = std::div(nidx - 4 - GetNodesPerEdge(mOrder) * 6, GetNodesPerFace(mOrder));
		IndexType face_gidx; int rotations; bool bucketflip;
		mConnectivityList[eidx].GetFaceNode(split.quot, &face_gidx, &rotations, &bucketflip);

		if (bucketflip)
		{
			split.rem = GetBucketflippedFaceIdx(split.rem, mOrder);
		}
		if (rotations > 0)
		{
			split.rem = GetRotatedFaceIdx(split.rem, mOrder, rotations);
		}
		return face_gidx + split.rem;
	}

	template<class IndexType>
	IndexType TetrahedralMesh<IndexType>::GetGlobalIndexBody(int eidx, int nidx)
	{
		int rem = nidx - 4 - GetNodesPerEdge(mOrder) * 6 - GetNodesPerFace(mOrder) * 4;
		return mConnectivityList[eidx].GetBodyNode() + rem;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetNumNodes() const
	{
		return this->mNodesInMesh;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetNumElements() const
	{
		return this->mElementsInMesh;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetNumNodesPerElement() const
	{
		return this->mNodesPrElement;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetNumLinearNodes()
	{
		return this->mOrgNoPoints;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetOrder() const
	{
		return this->mOrder;
	}

	template<class IndexType>
	typename TetrahedralMesh<IndexType>::IJKLType* TetrahedralMesh<IndexType>::GetIJKL(int lidx)
	{
		ASSERT(lidx >= 0);
		ASSERT(lidx < mNodesPrElement);

		return &mIJKL[MATRIX_INDEX(lidx, 0, 4)];
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetGlobalIndicesForSurfaceTriangle(const IndexType c1, const IndexType c2, const IndexType c3, IndexType *out, int* IJKs)
	{
		// Note that given that c1, c2, c3 are corner nodes of a surface, then it should be safe to assume that at most 1 element has all three nodes as corner nodes

		std::tuple<int, int, int> face_pairs[]{
			std::make_tuple(0, 1, 2), std::make_tuple(0, 1, 3),
			std::make_tuple(0, 2, 3), std::make_tuple(1, 2, 3)
		};

		int eidx = -1;
		int faceidx = -1;
		for (int i = 0; i < GetNumElements() && eidx < 0; i++)
		{
			CompressedElement element = mConnectivityList[i];
			for (int j = 0; j < 4; j++)
			{
				int f1, f2, f3; std::tie(f1, f2, f3) = face_pairs[j];
				if (c1 == element.GetCornerNode(f1) && c2 == element.GetCornerNode(f2) && c3 == element.GetCornerNode(f3) ||
					c1 == element.GetCornerNode(f1) && c3 == element.GetCornerNode(f2) && c2 == element.GetCornerNode(f3) ||
					c2 == element.GetCornerNode(f1) && c1 == element.GetCornerNode(f2) && c3 == element.GetCornerNode(f3) ||
					c2 == element.GetCornerNode(f1) && c3 == element.GetCornerNode(f2) && c1 == element.GetCornerNode(f3) ||
					c3 == element.GetCornerNode(f1) && c1 == element.GetCornerNode(f2) && c2 == element.GetCornerNode(f3) ||
					c3 == element.GetCornerNode(f1) && c2 == element.GetCornerNode(f2) && c1 == element.GetCornerNode(f3))
				{
					eidx = i;
					faceidx = j;
					break;
				}
			}
		}

		ASSERT(eidx >= 0);
		ASSERT(faceidx >= 0);

		int f1, f2, f3; std::tie(f1, f2, f3) = face_pairs[faceidx];

		out[0] = GetGlobalIndex(eidx, f1);
		// Only the 0th entry should be anything but zero, it should be equal to mOrder
		IJKs[0] = mIJKL[f1 * 4 + f1]; IJKs[1] = mIJKL[f1 * 4 + f2]; IJKs[2] = mIJKL[f1 * 4 + f3];
		ASSERT(IJKs[0] == mOrder);
		ASSERT(IJKs[1] == 0);
		ASSERT(IJKs[2] == 0);

		out[1] = GetGlobalIndex(eidx, f2);
		// Only the 1st entry should be anything but zero, it should be equal to mOrder
		IJKs[3] = mIJKL[f2 * 4 + f1]; IJKs[4] = mIJKL[f2 * 4 + f2]; IJKs[5] = mIJKL[f2 * 4 + f3];
		ASSERT(IJKs[3] == 0);
		ASSERT(IJKs[4] == mOrder);
		ASSERT(IJKs[5] == 0);

		out[2] = GetGlobalIndex(eidx, f3);
		// Only the 2nd entry should be anything but zero, it should be equal to mOrder
		IJKs[6] = mIJKL[f3 * 4 + f1]; IJKs[7] = mIJKL[f3 * 4 + f2]; IJKs[8] = mIJKL[f3 * 4 + f3];
		ASSERT(IJKs[6] == 0);
		ASSERT(IJKs[7] == 0);
		ASSERT(IJKs[8] == mOrder);

		int nextidx = 3;
		if (mOrder > 1)
		{
			// add the 3 edges
			std::tuple<int, int> edge_pairs[]{
				std::make_tuple(0, 1),
				std::make_tuple(0, 2),
				std::make_tuple(0, 3),
				std::make_tuple(1, 2),
				std::make_tuple(1, 3),
				std::make_tuple(2, 3)
			};

			// Edge 1
			int e1 = f1, e2 = f2;
			int edgeidx = -1;
			for (int i = 0; i < 6; i++)
			{
				std::tuple<int, int> edge = edge_pairs[i];
				if (e1 == std::get<0>(edge) && e2 == std::get<1>(edge) ||
					e2 == std::get<0>(edge) && e1 == std::get<1>(edge))
				{
					edgeidx = i;
					break;
				}
			}
			ASSERT(edgeidx >= 0);

			int nextEdgeIdx = 4 + GetNodesPerEdge(mOrder) * edgeidx;
			for (int i = 0; i < GetNodesPerEdge(mOrder); i++)
			{
				out[nextidx] = GetGlobalIndex(eidx, nextEdgeIdx);
				IJKs[nextidx * 3 + 0] = mIJKL[nextEdgeIdx * 4 + f1];
				IJKs[nextidx * 3 + 1] = mIJKL[nextEdgeIdx * 4 + f2];
				IJKs[nextidx * 3 + 2] = mIJKL[nextEdgeIdx * 4 + f3];

				nextidx++;
				nextEdgeIdx++;
			}

			// Edge 2
			e1 = f2, e2 = f3;
			edgeidx = -1;
			for (int i = 0; i < 6; i++)
			{
				std::tuple<int, int> edge = edge_pairs[i];
				if (e1 == std::get<0>(edge) && e2 == std::get<1>(edge) ||
					e2 == std::get<0>(edge) && e1 == std::get<1>(edge))
				{
					edgeidx = i;
					break;
				}
			}
			ASSERT(edgeidx >= 0);

			nextEdgeIdx = 4 + GetNodesPerEdge(mOrder) * edgeidx;
			for (int i = 0; i < GetNodesPerEdge(mOrder); i++)
			{
				out[nextidx] = GetGlobalIndex(eidx, nextEdgeIdx);
				IJKs[nextidx * 3 + 0] = mIJKL[nextEdgeIdx * 4 + f1];
				IJKs[nextidx * 3 + 1] = mIJKL[nextEdgeIdx * 4 + f2];
				IJKs[nextidx * 3 + 2] = mIJKL[nextEdgeIdx * 4 + f3];

				nextidx++;
				nextEdgeIdx++;
			}

			// Edge 3
			e1 = f3, e2 = f1;
			edgeidx = -1;
			for (int i = 0; i < 6; i++)
			{
				std::tuple<int, int> edge = edge_pairs[i];
				if (e1 == std::get<0>(edge) && e2 == std::get<1>(edge) ||
					e2 == std::get<0>(edge) && e1 == std::get<1>(edge))
				{
					edgeidx = i;
					break;
				}
			}
			ASSERT(edgeidx >= 0);

			nextEdgeIdx = 4 + GetNodesPerEdge(mOrder) * edgeidx;
			for (int i = 0; i < GetNodesPerEdge(mOrder); i++)
			{
				out[nextidx] = GetGlobalIndex(eidx, nextEdgeIdx);
				IJKs[nextidx * 3 + 0] = mIJKL[nextEdgeIdx * 4 + f1];
				IJKs[nextidx * 3 + 1] = mIJKL[nextEdgeIdx * 4 + f2];
				IJKs[nextidx * 3 + 2] = mIJKL[nextEdgeIdx * 4 + f3];

				nextidx++;
				nextEdgeIdx++;
			}
		}

		if (mOrder > 2)
		{
			// add the face

			int faceLocalidx = 4 + GetNodesPerEdge(mOrder) * 6 + GetNodesPerFace(mOrder) * faceidx;
			for (int i = 0; i < GetNodesPerFace(mOrder); i++)
			{
				out[nextidx] = GetGlobalIndex(eidx, faceLocalidx);
				IJKs[nextidx * 3 + 0] = mIJKL[faceLocalidx * 4 + f1];
				IJKs[nextidx * 3 + 1] = mIJKL[faceLocalidx * 4 + f2];
				IJKs[nextidx * 3 + 2] = mIJKL[faceLocalidx * 4 + f3];

				nextidx++;
				faceLocalidx++;
			}
		}
		return eidx;
	}

	// TODO the IJKL matrix code and build is exactly the same as the naive one, middleman?
	template<class IndexType>
	void TetrahedralMesh<IndexType>::BuildHigherOrderMesh(int noPoints, StridedVector<uint32> connectivity)
	{
		int nextIJKLRow = 0;
		IndexType nextGlobalIdx = noPoints;

		// -- Build the corner nodes --
		if (mOrder > 0)
		{
			// The 4 IJKL entries
			mIJKL[MATRIX_INDEX(nextIJKLRow++, 0, 4)] = mOrder;
			mIJKL[MATRIX_INDEX(nextIJKLRow++, 1, 4)] = mOrder;
			mIJKL[MATRIX_INDEX(nextIJKLRow++, 2, 4)] = mOrder;
			mIJKL[MATRIX_INDEX(nextIJKLRow++, 3, 4)] = mOrder;

			// The ConnectivityList entries
			for (int eidx = 0; eidx < mElementsInMesh; eidx++)
			{
				uint32 *gidxs = connectivity.GetAt(eidx);
				for (int i = 0; i < 4; i++)
				{
					IndexType gidx = (IndexType)(*(gidxs + i));
					mConnectivityList[eidx].SetCornerNode(i, gidx);

					auto it = std::find_if(mInterpolationscheme.begin(), mInterpolationscheme.end(),
						[gidx](const InterpolationTuple& e) { return std::get<0>(e) == gidx; });
					if (it == mInterpolationscheme.end())
					{
						mInterpolationscheme.push_back(std::make_tuple(gidx, eidx, i));
					}
				}
			}
		}

		// -- Build the edge nodes --
		if (mOrder > 1)
		{
			std::tuple<int, int> edge_pairs[]{
				std::make_tuple(0, 1),
				std::make_tuple(0, 2),
				std::make_tuple(0, 3),
				std::make_tuple(1, 2),
				std::make_tuple(1, 3),
				std::make_tuple(2, 3)
			};

			int no_nodes_pr_edge = GetNodesPerEdge(mOrder);


			// The IJKL entries
			for (int i = 0; i < 6; i++)
			{
				int e1, e2; std::tie(e1, e2) = edge_pairs[i];
				for (IJKLType k = mOrder - 1; k > 0; k--)
				{
					mIJKL[MATRIX_INDEX(nextIJKLRow, e1, 4)] = k;
					mIJKL[MATRIX_INDEX(nextIJKLRow++, e2, 4)] = mOrder - k;
				}
			}

			// The ConnectivityList entries

			// tuple<gidx of c1, gidx of c2, gidx of e1 of edge (c1,c2)>
			std::vector<tuple3IndexType> processed_edges;

			int lidx1, lidx2; IndexType gidx1, gidx2;
			for (int eidx = 0; eidx < mElementsInMesh; eidx++)
			{
				uint32 *gidxs = connectivity.GetAt(eidx);
				for (int i = 0; i < 6; i++)
				{
					std::tie(lidx1, lidx2) = edge_pairs[i];
					gidx1 = (IndexType)(*(gidxs + lidx1));
					gidx2 = (IndexType)(*(gidxs + lidx2));

					// First check if the edge has already been processed
					auto it = std::find_if(processed_edges.begin(), processed_edges.end(),
						[gidx1, gidx2](const tuple3IndexType& e) { return std::get<0>(e) == gidx1 && std::get<1>(e) == gidx2; });
					if (it != processed_edges.end())
					{
						IndexType gidx = std::get<2>(*it);
						mConnectivityList[eidx].SetEdgeNode(i, gidx, false);
						continue;
					}
					it = std::find_if(processed_edges.begin(), processed_edges.end(),
						[gidx1, gidx2](const tuple3IndexType& e) { return std::get<0>(e) == gidx2 && std::get<1>(e) == gidx1; });
					if (it != processed_edges.end())
					{
						IndexType gidx = std::get<2>(*it);
						mConnectivityList[eidx].SetEdgeNode(i, gidx + (no_nodes_pr_edge - 1), true);
						continue;
					}

					processed_edges.push_back(std::make_tuple(gidx1, gidx2, nextGlobalIdx));
					mConnectivityList[eidx].SetEdgeNode(i, nextGlobalIdx, false);

					int j_offset = 4 + i * no_nodes_pr_edge;
					for (int j = 0; j < no_nodes_pr_edge; j++)
					{
						mInterpolationscheme.push_back(std::make_tuple(nextGlobalIdx + j, eidx, j + j_offset));
					}
					nextGlobalIdx += no_nodes_pr_edge;
				}
			}
		}


		// -- Build the face nodes --
		if (mOrder > 2)
		{
			std::tuple<int, int, int> face_pairs[]{
				std::make_tuple(0, 1, 2), std::make_tuple(0, 1, 3),
				std::make_tuple(0, 2, 3), std::make_tuple(1, 2, 3)
			};

			int no_nodes_pr_face = GetNodesPerFace(mOrder);


			// The IJKL entries
			for (int i = 0; i < 4; i++)
			{
				int f1, f2, f3; std::tie(f1, f2, f3) = face_pairs[i];

				int minor_order = mOrder - 1;
				for (IJKLType v1 = minor_order; v1 > 0; v1--)
				{
					for (IJKLType v2 = minor_order - v1; v2 > 0; v2--)
					{
						IJKLType v3 = (mOrder - v1) - v2;

						mIJKL[MATRIX_INDEX(nextIJKLRow, f1, 4)] = v1;
						mIJKL[MATRIX_INDEX(nextIJKLRow, f2, 4)] = v2;
						mIJKL[MATRIX_INDEX(nextIJKLRow++, f3, 4)] = v3;
					}
				}
			}

			// The ConnectivityList entries
			// tuple<gidx of fidx1, gidx of fidx2, gidx of fidx3, gidx of f1 of face (fidx1,fidx2,fidx3)>
			std::vector<tuple4IndexType> processed_faces;

			auto checkProcessedFace = [&processed_faces, no_nodes_pr_face, this]
			(int eidx, int fidx, IndexType gidx1, IndexType gidx2, IndexType gidx3, int rotations, bool bucketflip)->bool
			{
				auto it = std::find_if(processed_faces.begin(), processed_faces.end(),
					[gidx1, gidx2, gidx3](const tuple4IndexType& e) {
					return std::get<0>(e) == gidx1 && std::get<1>(e) == gidx2 && std::get<2>(e) == gidx3;
				});

				if (it != processed_faces.end())
				{
					IndexType gidx = std::get<3>(*it);
					mConnectivityList[eidx].SetFaceNode(fidx, gidx, rotations, bucketflip);
					return true;
				}
				return false;
			};

			int lidx1, lidx2, lidx3; IndexType gidx1, gidx2, gidx3;
			for (int eidx = 0; eidx < mElementsInMesh; eidx++)
			{
				uint32 *gidxs = connectivity.GetAt(eidx);
				for (int i = 0; i < 4; i++)
				{
					std::tie(lidx1, lidx2, lidx3) = face_pairs[i];
					gidx1 = (IndexType)(*(gidxs + lidx1));
					gidx2 = (IndexType)(*(gidxs + lidx2));
					gidx3 = (IndexType)(*(gidxs + lidx3));

					// First check if the face has already been processed
					if (!(checkProcessedFace(eidx, i, gidx1, gidx2, gidx3, 0, false)
						|| checkProcessedFace(eidx, i, gidx3, gidx1, gidx2, 1, false)
						|| checkProcessedFace(eidx, i, gidx2, gidx3, gidx1, 2, false)
						|| checkProcessedFace(eidx, i, gidx1, gidx3, gidx2, 0, true)
						|| checkProcessedFace(eidx, i, gidx2, gidx1, gidx3, 1, true)
						|| checkProcessedFace(eidx, i, gidx3, gidx2, gidx1, 2, true)))
					{
						processed_faces.push_back(std::make_tuple(gidx1, gidx2, gidx3, nextGlobalIdx));
						mConnectivityList[eidx].SetFaceNode(i, nextGlobalIdx, 0, false);

						int j_offset = 4 + 6 * GetNodesPerEdge(mOrder) + i * no_nodes_pr_face;
						for (int j = 0; j < no_nodes_pr_face; j++)
						{
							mInterpolationscheme.push_back(std::make_tuple(nextGlobalIdx + j, eidx, j + j_offset));
						}
						nextGlobalIdx += no_nodes_pr_face;
					}
				}
			}
		}


		//-- Build the body nodes --
		if (mOrder > 3)
		{
			// The IJKL entries
			int body_order = mOrder - 4;

			for (IJKLType i = body_order + 1; i-- > 0;)
			{
				for (IJKLType j = std::max<IJKLType>(body_order - i, 0) + 1; j-- > 0;)
				{
					for (IJKLType k = std::max<IJKLType>(body_order - i - j, 0) + 1; k-- > 0;)
					{
						IJKLType l = body_order - i - j - k;
						mIJKL[MATRIX_INDEX(nextIJKLRow, 0, 4)] = i + 1;
						mIJKL[MATRIX_INDEX(nextIJKLRow, 1, 4)] = j + 1;
						mIJKL[MATRIX_INDEX(nextIJKLRow, 2, 4)] = k + 1;
						mIJKL[MATRIX_INDEX(nextIJKLRow++, 3, 4)] = l + 1;
					}
				}
			}

			// The ConnectivityList entries
			int no_nodes_pr_body = GetNodesPerBody(mOrder);

			int j_offset = 4 + 6 * GetNodesPerEdge(mOrder) + 4 * GetNodesPerFace(mOrder);
			for (int eidx = 0; eidx < mElementsInMesh; eidx++)
			{
				mConnectivityList[eidx].SetBodyNode(nextGlobalIdx);
				for (int j = 0; j < no_nodes_pr_body; j++)
				{
					mInterpolationscheme.push_back(std::make_tuple(nextGlobalIdx + j, eidx, j + j_offset));
				}
				nextGlobalIdx += no_nodes_pr_body;
			}
		}


		//// -- Debug printing --

		//// Debug-print the IJKL matrix
		//Printf("[I,J,K,L]\n");
		////Printf((std::to_string(order)).c_str());
		//for (int i = 0; i < mNodesPrElement; i++)
		//{
		//	Printf("[");
		//	Printf((std::to_string(mIJKL[MATRIX_INDEX(i, 0, 4)])).c_str());
		//	Printf(",");
		//	Printf((std::to_string(mIJKL[MATRIX_INDEX(i, 1, 4)])).c_str());
		//	Printf(",");
		//	Printf((std::to_string(mIJKL[MATRIX_INDEX(i, 2, 4)])).c_str());
		//	Printf(",");
		//	Printf((std::to_string(mIJKL[MATRIX_INDEX(i, 3, 4)])).c_str());
		//	Printf("]\n");
		//}
		//// -

		// //-

		//// Debug-print the ConnectivityList
		//Printf("ConnectivityList:\n");
		//for (int i = 0; i < mElementsInMesh; i++)
		//{
		//	CompressedElement ele = mConnectivityList[i];
		//	Printf("[");
		//	Printf((std::to_string(ele.c1)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.c2)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.c3)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.c4)).c_str()); Printf(",");

		//	Printf((std::to_string(ele.e1)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e1_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e2)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e2_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e3)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e3_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e4)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e4_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e5)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e5_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e6)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.e6_code)).c_str()); Printf(",");

		//	Printf((std::to_string(ele.f1)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f1_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f2)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f2_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f3)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f3_code)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f4)).c_str()); Printf(",");
		//	Printf((std::to_string(ele.f4_code)).c_str()); Printf(",");

		//	Printf((std::to_string(ele.b1)).c_str());
		//	Printf("]\n");
		//}

		//// Debug print the interpolation scheme
		//for (std::tuple<int, int, int> t : mInterpolationscheme)
		//{
		//	Printf("[gidx:");
		//	Printf((std::to_string(std::get<0>(t))).c_str());
		//	Printf(", eidx:");
		//	Printf((std::to_string(std::get<1>(t))).c_str());
		//	Printf(", lidx:");
		//	Printf((std::to_string(std::get<2>(t))).c_str());
		//	Printf("]\n");
		//}
	}

	// idx is expected to be the local idx of the face indices, not the actual global idx
	// ie idx \in [0, NumNodesPerFace(order)]
	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetRotatedFaceIdx(int idx, int order, int rotations)
	{
		int rotatedIdx = idx;
		int no_buckets = order - 2;

		// TODO, fix
		if (rotations == 2)
			rotations = 1;
		else if (rotations == 1)
			rotations = 2;

		while (rotations-- > 0)
		{
			int sum = 0;
			for (int i = 1; i <= no_buckets; i++)
			{
				sum += +i;
				if (rotatedIdx < sum)
				{
					rotatedIdx = Sum(1, no_buckets - (sum - rotatedIdx)) + no_buckets - i;
					break;
				}
			}
		}
		return rotatedIdx;
	}

	// idx is expected to be the local idx of the face indices, not the actual global idx
	// ie idx \in [0, NumNodesPerFace(order)]
	template<class IndexType>
	int TetrahedralMesh<IndexType>::GetBucketflippedFaceIdx(int idx, int order)
	{
		int no_buckets = order - 2;
		int sum = 0;
		for (int i = 1; i <= no_buckets; i++)
		{
			int next_sum = sum + i;
			if (idx < next_sum)
			{
				return sum + (next_sum - 1) - idx;
			}
			sum = next_sum;
		}
		return -1;
	}

	template<class IndexType>
	int TetrahedralMesh<IndexType>::Sum(int from, int to)
	{
		int sum = 0;
		for (int i = from; i <= to; i++)
		{
			sum += i;
		}
		return sum;
	}

} // namespace FEM_SYSTEM

#endif // TETRAHEDRAL_MESH_H
