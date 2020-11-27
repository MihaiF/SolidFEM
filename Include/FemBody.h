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

#pragma once

#include "FemDataStructures.h"
#include "MeshTools.h"

namespace FEM_SYSTEM
{
	class FemBody
	{
	public:
		void LoadFromXml(const char* path);
		void Prepare();
		
		void BuildBoundaryMesh();
		void UpdateBoundaryMesh(class FemPhysicsBase* femPhysics);
		void SaveBoundaryMesh(const char* path);
		
		bool HasVisualMesh() const { return !mVisualMesh.vertices.empty(); }
		void LoadVisualMesh(const char* path, const Vector3& offset, float scale);
		void MapVisualMesh();
		void UpdateVisualMesh();
		void SaveVisualMesh(const char* path);

		// non-const getters
		FemConfig& GetConfig() { return mConfig; }
		std::vector<Node>& GetNodes() { return mNodes; }
		std::vector<Tet>& GetTets() { return mTets; }
		std::vector<int>& GetFixedNodes() { return mFixedNodes; }

		// const getters
		const std::vector<Node>& GetNodes() const { return mNodes; }
		const std::vector<Tet>& GetTets() const { return mTets; }
		const Mesh& GetBoundaryMesh() const { return mBoundaryMesh; }
		const std::vector<int>& GetBoundaryNodes() const { return mBoundaryToNodesMap; }

	private:
		// volumetric mesh
		std::vector<Node> mNodes;
		std::vector<Tet> mTets;
		
		// setup
		std::vector<int> mFixedNodes;
		FemConfig mConfig;

		// boundary surface mesh
		Mesh mBoundaryMesh;
		std::vector<int> mNodesToBoundaryMap;
		std::vector<int> mBoundaryToNodesMap;
		std::vector<int> mBoundaryToElementsMap;

		// visual surface mesh
		Mesh mVisualMesh;
		std::vector<MeshInterp> mMapData;
	};
}