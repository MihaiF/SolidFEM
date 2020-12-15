#include "FemBody.h"
#include "FemPhysicsBase.h"
#include "FemIO.h"
#include "MeshTools.h"

namespace FEM_SYSTEM
{
	inline real ComputeTetVolume(const Vector3R& x0, const Vector3R& x1, const Vector3R& x2, const Vector3R& x3)
	{
		Vector3R d1 = x1 - x0;
		Vector3R d2 = x2 - x0;
		Vector3R d3 = x3 - x0;
		Matrix3R mat(d1, d2, d3);
		return mat.Determinant() / 6;
	}

	inline void ComputeBarycentric(const Vector3R& p, const Vector3R& x0, const Vector3R& x1, const Vector3R& x2, const Vector3R& x3,
		real& w0, real& w1, real& w2, real& w3)
	{
		real vol = ComputeTetVolume(x0, x1, x2, x3);
		real vol0 = ComputeTetVolume(p, x1, x2, x3);
		real vol1 = ComputeTetVolume(x0, p, x2, x3);
		real vol2 = ComputeTetVolume(x0, x1, p, x3);
		real vol3 = ComputeTetVolume(x0, x1, x2, p);
		w0 = vol0 / vol;
		w1 = vol1 / vol;
		w2 = vol2 / vol;
		w3 = vol3 / vol;
	}

	void CreateCable(const CableDescriptor& descriptor,
		const std::vector<Node>& nodes,
		const std::vector<Tet>& tets,
		FemPhysicsBase* femPhysics,
		real scale/* = 1*/)
	{
		Cable cable;
		int numSprings = descriptor.divs;
		real div = descriptor.length / numSprings;
		cable.mCableNodes.resize(numSprings + 1);
		cable.mCablePositions.resize(numSprings + 1);
		// generate a sequence of points along x axis
		for (int i = 0; i <= numSprings; i++)
		{
			Vector3R pos = i * div * descriptor.dir;
			pos += descriptor.offset;
			// register the point to the tet mesh
			real minScore = 1e20;
			int elem = -1;
			real minScorePos = 1e20;
			int elemPos = -1;
			for (uint32 e = 0; e < tets.size(); e++)
			{
				real score = 0;
				const Tet& tet = tets[e];

				Vector3R x0 = scale * nodes[tet.idx[0]].pos;
				Vector3R x1 = scale * nodes[tet.idx[1]].pos;
				Vector3R x2 = scale * nodes[tet.idx[2]].pos;
				Vector3R x3 = scale * nodes[tet.idx[3]].pos;
				real w0, w1, w2, w3;
				ComputeBarycentric(pos, x0, x1, x2, x3, w0, w1, w2, w3);
				real eps = -0.001;
				bool positive = w0 > eps && w1 > eps && w2 > eps && w3 > eps;

				for (int j = 0; j < 4; j++)
				{
					Vector3R node = 100.0 * nodes[tet.idx[j]].pos;
					score += (node - pos).LengthSquared();
				}
				if (score < minScore)
				{
					minScore = score;
					elem = e;
				}
				if (positive && score < minScorePos)
				{
					minScorePos = score;
					elemPos = e;
				}
			}
			if (elemPos >= 0)
				elem = elemPos;
			const Tet& tet = tets[elem];
			Vector3R x0 = scale * nodes[tet.idx[0]].pos;
			Vector3R x1 = scale * nodes[tet.idx[1]].pos;
			Vector3R x2 = scale * nodes[tet.idx[2]].pos;
			Vector3R x3 = scale * nodes[tet.idx[3]].pos;
			real w0, w1, w2, w3;
			ComputeBarycentric(pos, x0, x1, x2, x3, w0, w1, w2, w3);
			//Printf("%g, %g, %g, %g\n", w0, w1, w2, w3);
			real eps = 0.1;
			if (descriptor.useFreeCable && (w0 < -eps || w1 < -eps || w2 < -eps || w3 < -eps))
				cable.mCableNodes[i].elem = -elem;
			else
				cable.mCableNodes[i].elem = elem;
			cable.mCableNodes[i].bary.Set(w0, w1, w2);
			cable.mCablePositions[i] = (1 / scale) * pos;
		}

		cable.mActuation = 1;
		cable.mCableRestLength = div / scale;
		cable.mCableStiffness = descriptor.stiffness;

		femPhysics->AddCable(cable);
	}

	bool FemBody::LoadFromXml(const char* path)
	{
		std::vector<uint32> surfTris;
		std::string visualPath;
		float scale;		
		if (!IO::LoadFromXmlFile(std::string(path).c_str(), mNodes, mTets, mFixedNodes, surfTris, mConfig, scale, visualPath, mCables))
		{
			Printf("Failed to load file\n");
			return false;
		}
		Prepare();
		BuildBoundaryMesh();
		if (!visualPath.empty())
		{
			LoadVisualMesh(visualPath.c_str(), Vector3(), scale);
		}
		return true;
	}

	void FemBody::Prepare()
	{
		for (size_t i = 0; i < mNodes.size(); i++)
		{
			mNodes[i].pos0 = mNodes[i].pos;
		}
		for (size_t i = 0; i < mTets.size(); i++)
		{
			Tet& tet = mTets[i];
			const Vector3R& x0 = mNodes[tet.idx[0]].pos0;
			const Vector3R& x1 = mNodes[tet.idx[1]].pos0;
			const Vector3R& x2 = mNodes[tet.idx[2]].pos0;
			const Vector3R& x3 = mNodes[tet.idx[3]].pos0;
			Vector3R d1 = x1 - x0;
			Vector3R d2 = x2 - x0;
			Vector3R d3 = x3 - x0;
			Matrix3R mat(d1, d2, d3); // this is the reference shape matrix Dm [Sifakis][Teran03]
			real vol = (mat.Determinant()) / 6.f; // signed volume of the tet
			if (vol < 0)
			{
				// swap the first two indices so that the volume is positive next time
				std::swap(tet.idx[0], tet.idx[1]);
			}
		}
		for (uint32 i = 0; i < mFixedNodes.size(); i++)
		{
			int idx = mFixedNodes[i];
			mNodes[idx].invMass = 0.f;
		}
	}

	void FemBody::LoadVisualMesh(const char* path, const Vector3& offset, float scale)
	{
		ImportMeshFromOBJ(path, mVisualMesh, scale);
		MapVisualMesh();
	}

	void FemBody::MapVisualMesh()
	{
		mMapData.resize(mVisualMesh.vertices.size());
		for (uint j = 0; j < mVisualMesh.vertices.size(); j++)
		{
			float minDistance = FLT_MAX;
			int index = 0;
			const Vector3& vert = mVisualMesh.vertices[j];
			Vector3 closestP;
			int id1, id2, id3;
			float u, v, w;
			for (uint i = 0; i < mBoundaryMesh.GetNumTriangles(); i++)
			{
				int id1_t, id2_t, id3_t;
				id1_t = mBoundaryMesh.indices[i * 3];
				id2_t = mBoundaryMesh.indices[i * 3 + 1];
				id3_t = mBoundaryMesh.indices[i * 3 + 2];
				const Vector3& x1 = mBoundaryMesh.vertices[id1_t];
				const Vector3& x2 = mBoundaryMesh.vertices[id2_t];
				const Vector3& x3 = mBoundaryMesh.vertices[id3_t];
				Vector3 closestP_temp;
				closestP_temp = ClosestPtPointTriangle(vert, x1, x2, x3);
				Vector3 delta = vert - closestP_temp;
				float d = delta.Length();
				if (d < minDistance)
				{
					minDistance = d;
					index = i; // triangle index
					closestP = closestP_temp;
					id1 = id1_t;
					id2 = id2_t;
					id3 = id3_t;
				}
			}
			MeshInterp m;
			//m.distance = minDistance;
			m.id1 = id1;
			m.id2 = id2;
			m.id3 = id3;
			m.point = closestP;
			m.tet_index = index; // triangle index

			Vector3 v1 = DoubleToFloat(mBoundaryMesh.vertices[id1]);
			Vector3 v2 = DoubleToFloat(mBoundaryMesh.vertices[id2]);
			Vector3 v3 = DoubleToFloat(mBoundaryMesh.vertices[id3]);

			Vector3 normal = cross(v2 - v1, v3 - v1);
			normal.Normalize();
			m.distance = normal.Dot(vert - closestP);

			BarycentricCoordinates(closestP, v1, v2, v3, u, v, w);
			m.u = u;
			m.v = v;
			m.w = w;
			m.vert_index = j;
			m.original = vert;
			mMapData[j] = m;
		}
	}

	void FemBody::UpdateVisualMesh()
	{
		for (uint i = 0; i < mMapData.size(); i++)
		{
			const MeshInterp& m = mMapData[i];
			const Vector3& p1 = mBoundaryMesh.vertices[m.id1];
			const Vector3& p2 = mBoundaryMesh.vertices[m.id2];
			const Vector3& p3 = mBoundaryMesh.vertices[m.id3];

			Vector3 normal = cross(p2 - p1, p3 - p1);
			normal.Normalize();

			Vector3 p = p1 * m.w + p2 * m.v + p3 * m.u;//barycentric
			p += normal * m.distance;

			mVisualMesh.vertices[m.vert_index] = p;
		}
		mVisualMesh.ComputeNormals(); // slow but only way for now
	}

	void FemBody::BuildBoundaryMesh()
	{
		std::vector<Face> faces;
		for (uint32 ii = 0; ii < mTets.size(); ii++)
		{
			uint32 i = mTets[ii].idx[0];
			uint32 j = mTets[ii].idx[1];
			uint32 k = mTets[ii].idx[2];
			uint32 l = mTets[ii].idx[3];

			faces.push_back(Face(j, k, l, ii, 0));
			faces.push_back(Face(i, k, l, ii, 1));
			faces.push_back(Face(i, j, l, ii, 2));
			faces.push_back(Face(i, j, k, ii, 3));
		}

		std::sort(faces.begin(), faces.end(), CompareFaces);

		std::vector<Face> boundary;
		std::vector<uint16> boundaryVertices;
		for (size_t i = 0; i < faces.size(); i++)
		{
			if (i < faces.size() - 1 &&
				faces[i].i == faces[i + 1].i &&
				faces[i].j == faces[i + 1].j &&
				faces[i].k == faces[i + 1].k)
			{
				// if the next triangle is the same then this is clearly not a boundary triangle
				i++;
				continue;
			}
			boundary.push_back(faces[i]);
			boundaryVertices.push_back(faces[i].i);
			boundaryVertices.push_back(faces[i].j);
			boundaryVertices.push_back(faces[i].k);
		}

		// identify unique boundary vertices
		std::sort(boundaryVertices.begin(), boundaryVertices.end());

		// create vertex buffer and retain mappings
		std::vector<Vector3> vertices;
		mNodesToBoundaryMap.resize(mNodes.size());
		std::fill(mNodesToBoundaryMap.begin(), mNodesToBoundaryMap.end(), -1);
		mBoundaryToNodesMap.clear();
		for (size_t i = 0; i < boundaryVertices.size(); i++)
		{
			if (i > 0 && boundaryVertices[i] == boundaryVertices[i - 1])
				continue;
			vertices.push_back(DoubleToFloat(mNodes[boundaryVertices[i]].pos));
			int idx = (int)vertices.size() - 1;
			mNodesToBoundaryMap[boundaryVertices[i]] = idx;
			mBoundaryToNodesMap.push_back(boundaryVertices[i]);
		}

		// generate mesh
		mBoundaryMesh.vertices = vertices;
		const std::vector<int>& map = mNodesToBoundaryMap; // map from boundary mesh vertex indices to node indices

		// add the triangles
		mBoundaryMesh.indices.clear();
		mBoundaryToElementsMap.resize(boundary.size());
		for (size_t idx = 0; idx < boundary.size(); idx++)
		{
			// use the 4th vertex in the tet to identify which way is inside (for triangle winding)
			int ii = boundary[idx].elem;
			uint16 i = mTets[ii].idx[0];
			uint16 j = mTets[ii].idx[1];
			uint16 k = mTets[ii].idx[2];
			uint16 l = mTets[ii].idx[3];

			mBoundaryToElementsMap[idx] = ii;
			if (boundary[idx].face == 0)
			{
				const Vector3 v1 = DoubleToFloat(mNodes[j].pos);
				const Vector3 v2 = DoubleToFloat(mNodes[k].pos);
				const Vector3 v3 = DoubleToFloat(mNodes[l].pos);
				const Vector3 v4 = DoubleToFloat(mNodes[i].pos);
				Vector3 n = cross(v2 - v1, v3 - v1);
				bool flip = n.Dot(v4 - v1) > 0;
				mBoundaryMesh.AddTriangle(map[j], map[k], map[l], flip);
			}
			if (boundary[idx].face == 1)
			{
				const Vector3 v1 = DoubleToFloat(mNodes[i].pos);
				const Vector3 v2 = DoubleToFloat(mNodes[k].pos);
				const Vector3 v3 = DoubleToFloat(mNodes[l].pos);
				const Vector3 v4 = DoubleToFloat(mNodes[j].pos);
				Vector3 n = cross(v2 - v1, v3 - v1);
				bool flip = n.Dot(v4 - v1) > 0;
				mBoundaryMesh.AddTriangle(map[i], map[k], map[l], flip);
			}
			if (boundary[idx].face == 2)
			{
				const Vector3 v1 = DoubleToFloat(mNodes[i].pos);
				const Vector3 v2 = DoubleToFloat(mNodes[j].pos);
				const Vector3 v3 = DoubleToFloat(mNodes[l].pos);
				const Vector3 v4 = DoubleToFloat(mNodes[k].pos);
				Vector3 n = cross(v2 - v1, v3 - v1);
				bool flip = n.Dot(v4 - v1) > 0;
				mBoundaryMesh.AddTriangle(map[i], map[j], map[l], flip);
			}
			if (boundary[idx].face == 3)
			{
				const Vector3 v1 = DoubleToFloat(mNodes[i].pos);
				const Vector3 v2 = DoubleToFloat(mNodes[j].pos);
				const Vector3 v3 = DoubleToFloat(mNodes[k].pos);
				const Vector3 v4 = DoubleToFloat(mNodes[l].pos);
				Vector3 n = cross(v2 - v1, v3 - v1);
				bool flip = n.Dot(v4 - v1) > 0;
				mBoundaryMesh.AddTriangle(map[i], map[j], map[k], flip);
			}
		}

		mBoundaryMesh.ComputeNormals();
	}

	void FemBody::UpdateBoundaryMesh(FemPhysicsBase* femPhysics)
	{
		std::vector<real> scalarField;
		for (uint32 i = 0; i < mNodes.size(); i++)
		{
			Vector3R v = femPhysics->GetDeformedPosition(i);
			if (mNodesToBoundaryMap[i] >= 0)
			{
				mBoundaryMesh.vertices[mNodesToBoundaryMap[i]] = DoubleToFloat(v);
			}
		}
		mBoundaryMesh.ComputeNormals();
	}

	void FemBody::SaveBoundaryMesh(const char* path)
	{
		FILE* file;
		fopen_s(&file, path, "wt");
		ExportMeshToOBJ(file, mBoundaryMesh);
		fclose(file);
	}

	void FemBody::SaveVisualMesh(const char* path)
	{
		FILE* file;
		fopen_s(&file, path, "wt");
		ExportMeshToOBJ(file, mVisualMesh);
		fclose(file);
	}
}