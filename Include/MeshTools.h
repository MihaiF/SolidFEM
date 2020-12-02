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

#include <vector>
#include <Math/Vector3.h>
#include <Engine/Types.h>

using Math::Vector3;

struct Mesh
{
	std::vector<Vector3> vertices;
	std::vector<Vector3> normals;
	std::vector<Vector3> colors;
	//std::vector<Vector2> uvs;
	std::vector<uint32> indices;

	size_t GetNumTriangles() const { return indices.size() / 3; }

	void AddTriangle(uint32 a, uint32 b, uint32 c, bool flip = false)
	{
		indices.push_back(a);
		if (!flip)
		{
			indices.push_back(b);
			indices.push_back(c);
		}
		else
		{
			indices.push_back(c);
			indices.push_back(b);
		}
	}

	void ComputeNormals(bool flip = false);
};

struct Face
{
	uint32 i, j, k;
	int elem;
	int face;

	Face(uint32 a, uint32 b, uint32 c, int e, int f) : elem(e), face(f)
	{
		i = std::min(a, std::min(b, c));
		k = std::max(a, std::max(b, c));
		j = a + b + c - i - k;
	}
};

inline bool CompareFaces(const Face& t1, const Face& t2)
{
	if (t1.i < t2.i)
		return true;
	if (t1.i == t2.i)
	{
		if (t1.j < t2.j)
			return true;
		if (t1.j == t2.j && t1.k < t2.k)
			return true;
	}
	return false;
}

bool LoadMesh(const char* path, Mesh& collMesh, const Vector3& offset, float scale = 1, bool flipYZ = false);

void BarycentricCoordinates(const Vector3& p,
	const Vector3& a,
	const Vector3& b,
	const Vector3& c,
	float& u,
	float& v,
	float& w);

Vector3 ClosestPtPointTriangle(const Vector3& p,
	const Vector3& a,
	const Vector3& b,
	const Vector3& c);

void ExportMeshToOBJ(FILE* txt, const Mesh& mesh);

void ImportMeshFromOBJ(const char* path, Mesh& mesh, float scale = 1);


