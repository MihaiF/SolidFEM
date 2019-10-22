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

#ifndef MESH_FACTORY_H
#define MESH_FACTORY_H
// TODO bad naming, 'MeshFactory'? As it also interpolates, and pot. other mesh-helper functions


// Custom references
#include "StridedVector.h"
#include <Engine/Types.h>
#include <Engine/Utils.h>

namespace FEM_SYSTEM
{
	class MeshFactory {
	public:
		template<class MeshType>
		static MeshType* MakeMesh(uint order, uint noPoints, StridedVector<uint16> connectivity, uint noElements)
		{
			// TODO assert that connectivity is not null, how to?

			MeshType* mesh = new MeshType();
			MakeMesh(mesh, order, noPoints, connectivity, noElements);
			return mesh;
		}

		template<class MeshType>
		static void MakeMesh(MeshType* mesh, uint order, uint noPoints, StridedVector<uint16> connectivity, uint noElements)
		{
			// TODO assert that connectivity is not null, how to?
			// TODO assert that mesh is not null, how to?

			mesh->InstanceWith(order, noPoints, connectivity, noElements);
		}


		// TODO move into some util file?
		template<class T, class IJKLType>
		static T LerpTetrahedra(T c1, T c2, T c3, T c4, IJKLType I, IJKLType J, IJKLType K, IJKLType L, int order)
		{
			ASSERT(order > 0);
			ASSERT(I >= 0);
			ASSERT(J >= 0);
			ASSERT(K >= 0);
			ASSERT(L >= 0);
			ASSERT((int)I <= order);
			ASSERT((int)J <= order);
			ASSERT((int)K <= order);
			ASSERT((int)L <= order);

			float forder = (float)order;
			return LerpTetrahedra(c1, c2, c3, c4, I / forder, J / forder, K / forder, L / forder);
		}

		// TODO move into some util file?
		template<class T>
		static T LerpTetrahedra(T c1, T c2, T c3, T c4, float L1, float L2, float L3, float L4)
		{
			return c1 * L1 + c2 * L2 + c3 * L3 + c4 * L4;
		}

		template<class MeshType, class DataType>
		static void Interpolate(MeshType &mesh, DataType* in, DataType* out)
		{
			// TODO assert that mesh is not null, how to?
			// ASSERT that in has the expected length, can you do that with a pointer?
			// ASSERT that out has the expected length, can you do that with a pointer?

			typename MeshType::mIndexType gidx; int eidx, lidx;
			typename MeshType::IJKLType I, J, K, L;
			DataType c1, c2, c3, c4;
			typename MeshType::IJKLType* IJKL;

			for (int i = 0; i < mesh.GetNumNodes(); i++)
			{
				std::tie(gidx, eidx, lidx) = mesh.mInterpolationscheme[i];

				IJKL = mesh.GetIJKL(lidx);
				I = IJKL[0];
				J = IJKL[1];
				K = IJKL[2];
				L = IJKL[3];

				c1 = in[mesh.GetGlobalIndex(eidx, 0)];
				c2 = in[mesh.GetGlobalIndex(eidx, 1)];
				c3 = in[mesh.GetGlobalIndex(eidx, 2)];
				c4 = in[mesh.GetGlobalIndex(eidx, 3)];

				out[gidx] = LerpTetrahedra(c1, c2, c3, c4, I, J, K, L, mesh.GetOrder());
			}
		}

		template<class MeshType, class DataType>
		static void Interpolate(MeshType &mesh, StridedVector<DataType> in, DataType* out)
		{
			// TODO assert that mesh is not null, how to?
			// ASSERT that in has the expected length, can you do that with a pointer?
			// ASSERT that out has the expected length, can you do that with a pointer?

			typename MeshType::mIndexType gidx; int eidx, lidx;
			typename MeshType::IJKLType I, J, K, L;
			DataType c1, c2, c3, c4;
			typename MeshType::IJKLType* IJKL;

			for (int i = 0; i < mesh.GetNumNodes(); i++)
			{
				std::tie(gidx, eidx, lidx) = mesh.mInterpolationscheme[i];

				IJKL = mesh.GetIJKL(lidx);
				I = IJKL[0];
				J = IJKL[1];
				K = IJKL[2];
				L = IJKL[3];

				c1 = *in.GetAt((uint)mesh.GetGlobalIndex(eidx, 0));
				c2 = *in.GetAt((uint)mesh.GetGlobalIndex(eidx, 1));
				c3 = *in.GetAt((uint)mesh.GetGlobalIndex(eidx, 2));
				c4 = *in.GetAt((uint)mesh.GetGlobalIndex(eidx, 3));

				out[gidx] = LerpTetrahedra(c1, c2, c3, c4, I, J, K, L, mesh.GetOrder());
			}
		}
	};
}

#endif // MESH_FACTORY_H
