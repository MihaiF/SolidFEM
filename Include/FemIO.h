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

#ifndef FEM_IO_H
#define FEM_IO_H

#include "FemDataStructures.h"
#include <set>
#include <fstream>

enum SwapCoordsType
{
	SCT_SWAP_NONE,
	SCT_SWAP_X_AND_Y,
	SCT_SWAP_X_AND_Z,
	SCT_SWAP_Y_AND_Z,
};

namespace FEM_SYSTEM
{
	class FemPhysicsBase;
	namespace IO
	{
		bool LoadFromFebFile(const char* path, std::vector<Node>& nodes, std::vector<Tet>& tets,
			std::vector<int>& fixedNodes, std::vector<uint32>& surfTris, std::set<uint32>& innerSurface, 
			real scale, FemConfig* config = nullptr, int* bcFlag = nullptr, std::vector<uint32>* bcIndices = nullptr);
		bool LoadFromXmlFile(const char* path, std::vector<Node>& nodes, std::vector<Tet>& tets,
			std::vector<int>& fixedNodes, std::vector<uint32>& surfTris, FemConfig& config);
		bool LoadFromVolFile(const char* path, real scale, std::vector<Node>& nodes, std::vector<Tet>& tets);
		bool LoadFromTet1File(const char* path, real scale, const Vector3R& offset, std::vector<Node>& nodes, std::vector<Tet>& tets);
		void ExportToVTKFile(std::fstream& outfile, FemPhysicsBase* femPhysics);
		void ExportToVTKHeatMap(std::fstream& outfile, FemPhysicsBase* femPhysics);
	}
}

#endif // FEM_IO_H
