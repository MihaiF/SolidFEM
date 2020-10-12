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

#ifndef LINEAR_TETRAHEDRAL_MESH_H
#define LINEAR_TETRAHEDRAL_MESH_H

// C++ references
#include <vector>

// Custom references
#include "StridedVector.h"
#include "FemDataStructures.h"
#include "ITetrahedralMesh.h"
#include <Engine/Types.h>

namespace FEM_SYSTEM
{
	class LinearTetrahedralMesh : public ITetrahedralMesh
	{
	public:
		typedef uint32 mIndexType;
		typedef int IJKLType;

		void InstanceWith(int order, int noPoints, StridedVector<uint32> connectivity, int noElements)
		{
			mNumNodes = noPoints;
			mNumElems = noElements;
			mTets.resize(noElements);
			for (int i = 0; i < noElements; i++)
			{
				uint32* gidxs = connectivity.GetAt(i);
				for (int j = 0; j < 4; j++)
					mTets[i].idx[j] = gidxs[j];
			}
		}
		int GetNumNodes() const { return mNumNodes; }
		int GetNumNodesPerElement() const { return 4; }
		int GetNumElements() const { return mNumElems; }
		int GetNodesPerEdge(int) const { return 0; }
		int GetNodesPerFace(int) const { return 0; }
		int GetGlobalIndex(int eidx, int lidx) { return mTets[eidx].idx[lidx]; }
		int* GetIJKL(int lidx) { return mIJKL[lidx]; }
		int GetOrder() const { return 1; }

	private:
		uint32 mNumNodes;
		uint32 mNumElems;
		std::vector<Tet> mTets;
		int mIJKL[4][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
	};
}

#endif // LINEAR_TETRAHEDRAL_MESH_H