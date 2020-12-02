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

#include "MeshTools.h"

#include <Engine/Utils.h>

#include <assimp/cimport.h>        // Plain-C interface
#include <assimp/scene.h>          // Output data structure
#include <assimp/postprocess.h>    // Post processing flags

#include <iostream>
#include <fstream>
#include <string>
#include <regex>

void BarycentricCoordinates(const Vector3& p,
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

Vector3 ClosestPtPointTriangle(const Vector3& p,
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
	float vc = d1 * d4 - d3 * d2;
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
	float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
		float w = d2 / (d2 - d6);
		return a + w * ac; // barycentric coordinates (1-w,0,w)
	}
	// Check if P in edge region of BC, if so return projection of P onto BC
	float va = d3 * d6 - d5 * d4;
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

void CreateCollisionMesh(const struct aiScene* sc, const struct aiNode* nd, Mesh& collMesh,
	const Vector3& offset, float scale, bool flipYZ)
{
	ASSERT(sc != NULL);

	for (size_t n = 0; n < nd->mNumMeshes; ++n)
	{
		// found one
		const struct aiMesh* mesh = sc->mMeshes[nd->mMeshes[n]];

		size_t baseIdx = collMesh.vertices.size();
		for (size_t i = 0; i < mesh->mNumVertices; i++)
		{
			const aiVector3D& src = mesh->mVertices[i];
			aiVector3D t = /*nd->mTransformation * */src;
			Vector3 dst(t.x * scale, t.y * scale, t.z * scale);
			if (flipYZ)
				dst.Set(t.x * scale, t.z * scale, t.y * scale);
			const aiVector3D& srcN = mesh->mNormals[i];
			Vector3 n(srcN.x, srcN.y, srcN.z);
			n.Normalize();
			collMesh.vertices.push_back(dst + offset);
			collMesh.normals.push_back(n);

			if (mesh->GetNumColorChannels() > 0 && mesh->HasVertexColors(0))
			{
				Vector3 col(mesh->mColors[0][i].r, mesh->mColors[0][i].g, mesh->mColors[0][i].b);
				collMesh.colors.push_back(col);
			}

			//if (mesh->GetNumUVChannels() > 0 && mesh->HasTextureCoords(0))
			//{
			//	Vector2 uv(mesh->mTextureCoords[0][i].x, mesh->mTextureCoords[0][i].y);
			//	collMesh.uvs.push_back(uv);
			//}
		}

		for (size_t t = 0; t < mesh->mNumFaces; ++t)
		{
			const struct aiFace* face = &mesh->mFaces[t];
			for (size_t i = 0; i < face->mNumIndices; i += 3)
			{
				uint32 index0 = (uint32)(face->mIndices[i] + baseIdx);
				uint32 index1 = (uint32)(face->mIndices[i + 1] + baseIdx);
				uint32 index2 = (uint32)(face->mIndices[i + 2] + baseIdx);

				//Need this to avoid back face culling issue.
				if (flipYZ) {
					collMesh.indices.push_back(index0);
					collMesh.indices.push_back(index2);
					collMesh.indices.push_back(index1);
				}
				else {
					collMesh.indices.push_back(index0);
					collMesh.indices.push_back(index1);
					collMesh.indices.push_back(index2);
				}
			}
		}

		//return;
	}

	// recurse if no mesh found
	for (size_t n = 0; n < nd->mNumChildren; ++n)
	{
		CreateCollisionMesh(sc, nd->mChildren[n], collMesh, offset, scale, flipYZ);
	}
}

bool LoadMesh(const char* path, Mesh& collMesh, const Vector3& offset, float scale, bool flipYZ)
{
	const aiScene* scene = aiImportFile(path, 0);
	if (scene == NULL)
		return false;
	collMesh.indices.clear();
	collMesh.vertices.clear();
	CreateCollisionMesh(scene, scene->mRootNode, collMesh, offset, scale, flipYZ);
	return true;
}

void ExportMeshToText(FILE* txt, const Mesh& mesh)
{
	for (size_t j = 0; j < mesh.vertices.size(); j++)
	{
		fprintf(txt, "%f, %f, %f,\n", mesh.vertices[j].X(), mesh.vertices[j].Y(), mesh.vertices[j].Z());
	}
	for (size_t j = 0; j < mesh.normals.size(); j++)
	{
		fprintf(txt, "%f, %f, %f,\n", mesh.normals[j].X(), mesh.normals[j].Y(), mesh.normals[j].Z());
	}
	for (size_t j = 0; j < mesh.indices.size(); j += 3)
	{
		fprintf(txt, "%d, %d, %d,\n", mesh.indices[j], mesh.indices[j + 1], mesh.indices[j + 2]);
	}
}

void ExportMeshToOBJ(FILE* txt, const Mesh& mesh)
{
	for (size_t j = 0; j < mesh.vertices.size(); j++)
	{
		fprintf(txt, "vn %f %f %f\n", mesh.normals[j].X(), mesh.normals[j].Y(), mesh.normals[j].Z());
		fprintf(txt, "v %f %f %f\n", mesh.vertices[j].X(), mesh.vertices[j].Y(), mesh.vertices[j].Z());
	}
	for (size_t j = 0; j < mesh.indices.size(); j += 3)
	{
		fprintf(txt, "f %d//%d %d//%d %d//%d\n", mesh.indices[j] + 1, mesh.indices[j] + 1, mesh.indices[j + 1] + 1, mesh.indices[j + 1] + 1, mesh.indices[j + 2] + 1, mesh.indices[j + 2] + 1);
	}
}

void ImportMeshFromOBJ(const char* path, Mesh& mesh, float scale)
{
	mesh.vertices.clear();
	mesh.normals.clear();
	mesh.indices.clear();

	std::ifstream file(path);
	std::string str;
	while (std::getline(file, str)) 
	{
		// Regex for tokenizing whitespaces
		std::regex reg("\\s+");
		// Get an iterator after filtering through the regex
		std::sregex_token_iterator iter(str.begin(), str.end(), reg, -1);
		// Keep a dummy end iterator - Needed to construct a vector  using (start, end) iterators.
		std::sregex_token_iterator end;
		std::vector<std::string> tokens(iter, end);
		//for (auto a : tokens)
		//{
		//	std::cout << a << std::endl;
		//}

		if (tokens[0].compare("v") == 0)
		{
			float x = std::stof(tokens[1]);
			float y = std::stof(tokens[2]);
			float z = std::stof(tokens[3]);
			mesh.vertices.push_back(scale * Vector3(x, y, z));
		}

		if (tokens[0].compare("vn") == 0)
		{
			float x = std::stof(tokens[1]);
			float y = std::stof(tokens[2]);
			float z = std::stof(tokens[3]);
			mesh.normals.push_back(Vector3(x, y, z));
		}
		
		if (tokens[0].compare("f") == 0)
		{
			int i = std::stoi(tokens[1]);
			int j = std::stoi(tokens[2]);
			int k = std::stoi(tokens[3]);
			mesh.indices.push_back(i - 1);
			mesh.indices.push_back(j - 1);
			mesh.indices.push_back(k - 1);
		}
	}
}

void Mesh::ComputeNormals(bool flip /*= false*/)
{
	normals.resize(vertices.size());
	for (size_t i = 0; i < vertices.size(); i++)
		normals[i].SetZero();
	for (size_t i = 0; i < indices.size(); i += 3)
	{
		const int i1 = indices[i];
		const int i2 = indices[i + 1];
		const int i3 = indices[i + 2];
		const Vector3& v1 = vertices[i1];
		const Vector3& v2 = vertices[i2];
		const Vector3& v3 = vertices[i3];
		Vector3 n = (v2 - v1).Cross(v3 - v1);
#ifdef ANGLE_WEIGHTED
		n.Normalize();
		Vector3 e1 = v2 - v1;
		Vector3 e2 = v3 - v2;
		Vector3 e3 = v1 - v3;
		e1.Normalize();
		e2.Normalize();
		e3.Normalize();
		float t1 = acosf(-(e1 * e3));
		float t2 = acosf(-(e1 * e2));
		float t3 = acosf(-(e2 * e3));
		normals[i1] += (t1 / PI) * n;
		normals[i2] += (t2 / PI) * n;
		normals[i3] += (t3 / PI) * n;
#else
		normals[i1] += n;
		normals[i2] += n;
		normals[i3] += n;
#endif
	}
	for (size_t i = 0; i < vertices.size(); i++)
	{
		normals[i].Normalize();
		if (flip)
			normals[i].Flip();
	}
}
