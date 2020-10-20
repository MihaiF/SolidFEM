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

#include "FemIO.h"
#include "FemPhysicsBase.h"
#include "ElasticEnergy.h"
#include <tinyxml2.h>

namespace FEM_SYSTEM
{
	bool IO::LoadFromFebFile(const char* path, std::vector<Node>& nodes, std::vector<Tet>& tets, 
		std::vector<int>& fixedNodes, std::vector<uint32>& surfTris, std::set<uint32>& innerSurface, 
		real scale, FemConfig* femConfig, int* bcFlag, std::vector<uint32>* bcIndices)
	{
		fixedNodes.clear();
		nodes.clear();
		tets.clear();
		surfTris.clear();
		innerSurface.clear();

		typedef std::vector<uint32> NodeSet;
		std::map<std::string, NodeSet> nodeSets;
		std::map<std::string, std::string> nodeSetAliases;

		tinyxml2::XMLDocument doc;		
		if (doc.LoadFile(path) == tinyxml2::XML_SUCCESS)
		{
			//char str[128];
			tinyxml2::XMLElement* root = doc.FirstChildElement("febio_spec");
			if (!root)
				return false;

			tinyxml2::XMLElement* xGeom = root->FirstChildElement("Geometry");
			if (xGeom)
			{
				tinyxml2::XMLElement* xPart = xGeom->FirstChildElement("Part");
				if (xPart)
				{
					const char* name = xPart->Attribute("name");
					Printf("Part: %s\n", name);

					// node set aliases
					tinyxml2::XMLElement* xNodeSet = xGeom->FirstChildElement("NodeSet");
					while (xNodeSet)
					{
						const char* xName = xNodeSet->Attribute("name");
						if (xName)
						{
							Printf("NodeSet alias: %s\n", xName);
							tinyxml2::XMLElement* xAlias = xNodeSet->FirstChildElement("NodeSet");
							if (xAlias)
							{
								const char* alias = xAlias->Attribute("node_set");
								nodeSetAliases.insert({ xName, alias });
							}
						}

						xNodeSet = xNodeSet->NextSiblingElement("NodeSet");
					}

					xGeom = xPart; // hack to make the geometry reading work
				}
				// nodes
				tinyxml2::XMLElement* xNodes = xGeom->FirstChildElement("Nodes");
				if (xNodes)
				{
					tinyxml2::XMLElement* xNode = xNodes->FirstChildElement();
					while (xNode)
					{						
						const char* text = xNode->GetText();
						//Printf("%s\n", text);

						std::string str(text);
						size_t next1, next2;
						double x1 = std::stod(str, &next1);
						double x2 = std::stod(str.substr(next1 + 1), &next2);
						double x3 = std::stod(str.substr(next1 + next2 + 2));
						//Printf("%g %g %g\n", x1, x2, x3);

						Node node;
						node.pos = Vector3R(x1, x2, x3);
						node.pos.Scale(scale);
						//node.pos0 = node.pos;
						nodes.push_back(node);

						xNode = xNode->NextSiblingElement();
					}
				}

				// tets
				tinyxml2::XMLElement* xElements = xGeom->FirstChildElement("Elements");
				if (xElements)
				{
					tinyxml2::XMLElement* xElement = xElements->FirstChildElement();
					while (xElement)
					{
						const char* text = xElement->GetText();
						//Printf("%s\n", text);

						std::string str(text);
						size_t next1, next2, next3;
						int i1 = std::stoi(str, &next1);
						int i2 = std::stoi(str.substr(next1 + 1), &next2);
						int i3 = std::stoi(str.substr(next1 + next2 + 2), &next3);
						int i4 = std::stoi(str.substr(next1 + next2 + next3 + 3));
						//Printf("%d %d %d %d\n", i1, i2, i3, i4);

						Tet tet;
						tet.idx[0] = i1 - 1;
						tet.idx[1] = i2 - 1;
						tet.idx[2] = i3 - 1;
						tet.idx[3] = i4 - 1;
						tets.push_back(tet);

						xElement = xElement->NextSiblingElement();
					}
				}

				// node sets
				tinyxml2::XMLElement* xNodeSet = xGeom->FirstChildElement("NodeSet");
				while (xNodeSet)
				{
					const char* xName = xNodeSet->Attribute("name");
					if (xName)
					{
						Printf("NodeSet: %s\n", xName);
					}

					NodeSet ns;
					tinyxml2::XMLElement* xNode = xNodeSet->FirstChildElement();
					while (xNode)
					{
						if (xNode)
						{
							const char* xId = xNode->Attribute("id");
							if (xId)
							{
								int id = atoi(xId);
								ns.push_back(id - 1);
								//Printf("%d,", id);
							}
						}
						xNode = xNode->NextSiblingElement();
					}
					//Printf("\n");
					nodeSets.insert({ xName, ns });

					xNodeSet = xNodeSet->NextSiblingElement("NodeSet");
				}

				// pressure load
				tinyxml2::XMLElement* xSurface = xGeom->FirstChildElement("Surface");
				if (xSurface)
				{
					const char* xName = xSurface->Attribute("name");
					if (xName)
					{
						Printf("Surface: %s\n", xName);
					}

					tinyxml2::XMLElement* xTri = xSurface->FirstChildElement();
					while(xTri)
					{
						const char* text = xTri->GetText();
						//Printf("%s\n", text);

						std::string str(text);
						size_t next1, next2;
						int i1 = std::stoi(str, &next1);
						int i2 = std::stoi(str.substr(next1 + 1), &next2);
						int i3 = std::stoi(str.substr(next1 + next2 + 2));
						//Printf("%d %d %d\n", i1, i2, i3);

						i1--;
						i2--;
						i3--;
						surfTris.push_back(i1);
						surfTris.push_back(i2);
						surfTris.push_back(i3);

						innerSurface.insert(i1);
						innerSurface.insert(i2);
						innerSurface.insert(i3);

						xTri = xTri->NextSiblingElement();
					}
				}

			}

			tinyxml2::XMLElement* xBoundary = root->FirstChildElement("Boundary");
			if (xBoundary)
			{
				tinyxml2::XMLElement* xFix = xBoundary->FirstChildElement("fix");
				while (xFix)
				{
					const char* bc = xFix->Attribute("bc");
					const char* name = xFix->Attribute("node_set");

					NodeSet nodeSet;
					if (nodeSetAliases.empty())
					{
						nodeSet = nodeSets[name];
					}
					else
					{
						// TODO: use find
						std::string alias = nodeSetAliases[name];
						size_t point = alias.find('.');
						alias = alias.substr(point + 1);
						nodeSet = nodeSets[alias];
					}
					if (strcmp(bc, "x,y,z") == 0)
					{
						fixedNodes.insert(fixedNodes.end(), nodeSet.begin(), nodeSet.end());
					}
					else
					{
						*bcFlag = FemPhysicsBase::AXIS_X | FemPhysicsBase::AXIS_Y;
						bcIndices->insert(bcIndices->end(), nodeSet.begin(), nodeSet.end());
					}
					
					
					xFix = xFix->NextSiblingElement("fix");				
				}
			}

			tinyxml2::XMLElement* xLoads = root->FirstChildElement("Loads");
			if (xLoads)
			{
				tinyxml2::XMLElement* xSurfLoad = xLoads->FirstChildElement("surface_load");
				if (xSurfLoad)
				{
					tinyxml2::XMLElement* xPressure = xSurfLoad->FirstChildElement("pressure");
					if (xPressure)
					{
						const char* text = xPressure->GetText();
						int p = atoi(text);
						//Printf("%d\n", p);
						if (femConfig)
							femConfig->mAppliedPressure = p;
					}
				}

				tinyxml2::XMLElement* xBodyLoad = xLoads->FirstChildElement("body_load");
				if (xBodyLoad)
				{
					tinyxml2::XMLElement* xY = xBodyLoad->FirstChildElement("y");
					real gy = atof(xY->GetText());
					tinyxml2::XMLElement* xZ = xBodyLoad->FirstChildElement("z");
					real gz = atof(xZ->GetText());
					if (femConfig)
						femConfig->mGravity = gy; // TODO: z up
				}
			}

			tinyxml2::XMLElement* xMaterial = root->FirstChildElement("Material");
			if (xMaterial)
			{
				tinyxml2::XMLElement* xMat = xMaterial->FirstChildElement("material");
				if (xMat)
				{
					tinyxml2::XMLElement* xYoung = xMat->FirstChildElement("E");
					if (xYoung)
					{
						const char* text = xYoung->GetText();
						double E = atof(text);
						if (femConfig)
							femConfig->mYoungsModulus = E;
					}
					tinyxml2::XMLElement* xPoisson = xMat->FirstChildElement("v");
					if (xPoisson)
					{
						const char* text = xPoisson->GetText();
						double v = atof(text);
						if (femConfig)
							femConfig->mPoissonRatio = v;
					}
				}
			}

			return true;
		}

		return false;
	}

	void IO::ExportToVTKHeatMap(std::fstream& outfile, FemPhysicsBase* femPhysics)
	{
		IO::ExportToVTKFile(outfile, femPhysics);

		// displacements
		outfile << "POINT_DATA " << femPhysics->GetNumNodes() << std::endl;
		outfile << "VECTORS displacement float" << std::endl;

		for (uint32 i = 0; i < femPhysics->GetNumNodes(); i++)
		{
			auto u = femPhysics->GetDeformedPosition(i) - femPhysics->GetInitialPosition(i);
			outfile << std::to_string(u.x);
			outfile << " ";
			outfile << std::to_string(u.y);
			outfile << " ";
			outfile << std::to_string(u.z);
			outfile << "\n";
		}
		outfile << "\n";

		// stress
		outfile << "CELL_DATA " << femPhysics->GetNumElements() << std::endl;
		outfile << "TENSORS stress float" << std::endl;
		for (uint32 i = 0; i < femPhysics->GetNumElements(); i++)
		{
			auto P = ElasticEnergy::ComputeElementStress(femPhysics, i);
			outfile << P(0, 0) << " " << P(0, 1) << " " << P(0, 2) << std::endl;
			outfile << P(1, 0) << " " << P(1, 1) << " " << P(1, 2) << std::endl;
			outfile << P(2, 0) << " " << P(2, 1) << " " << P(2, 2) << std::endl;
			outfile << std::endl;
		}
		outfile << std::endl;
	}

	void IO::ExportToVTKFile(std::fstream& outfile, FemPhysicsBase* femPhysics)
	{
		// TODO it currently only supports meshes stored as TetrahedralMesh and uses GetInitialPosition
		int prefix_index = 0;

		outfile << "# vtk DataFile Version 2.0\n";
		outfile << "Unstructured Grid\n";
		outfile << "ASCII\n";
		outfile << "DATASET UNSTRUCTURED_GRID\n";

		outfile << "POINTS " << std::to_string(femPhysics->GetNumNodes()) << " double\n";
		for (uint32 i = 0; i < femPhysics->GetNumNodes(); i++)
		{
			auto v = femPhysics->GetDeformedPosition(i);
			outfile << std::to_string(v.x);
			outfile << " ";
			outfile << std::to_string(v.y);
			outfile << " ";
			outfile << std::to_string(v.z);
			outfile << "\n";
		}
		outfile << "\n";

		uint32 n = femPhysics->GetNumLocalNodes();
		outfile << "CELLS " << femPhysics->GetNumElements() << " " << femPhysics->GetNumElements() * (n + 1) << "\n";
		const uint32 map[10] = { 0, 1, 2, 3, 4, 7, 5, 6, 8, 9 }; // map from our local IJKL order convention to VTK
		for (uint32 i = 0; i < femPhysics->GetNumElements(); i++)
		{
			outfile << std::to_string(n);
			outfile << " ";
			for (uint32 j = 0; j < n; j++)
			{
				outfile << std::to_string(femPhysics->GetGlobalOriginalIndex(i, map[j]) + prefix_index);
				outfile << " ";
			}
			outfile << "\n";
		}
		outfile << "\n";

		outfile << "CELL_TYPES " << femPhysics->GetNumElements() << "\n";
		int VTK_TETRA = 10;
		for (uint32 i = 0; i < femPhysics->GetNumElements(); i++)
		{
			outfile << std::to_string(n > 4 ? 78 : VTK_TETRA);
			outfile << "\n";
		}
		outfile << "\n";
	}

	bool IO::LoadFromVolFile(const char* path, real scale, std::vector<Node>& nodes, std::vector<Tet>& tets)
	{
		FILE* ptr_file;
		char buf[1000];

		const int STATE_READ = 0;
		const int STATE_VOL = 1;
		const int STATE_ELEM = 2;
		const int STATE_POINTS = 3;
		const int STATE_POINT = 4;

		errno_t ret = fopen_s(&ptr_file, path, "r");
		if (ret != 0)
		{
			Printf("Couldn't open file %s\n", path);
			return false;
		}

		int state = STATE_READ;
		int nod = 0;
		while (fgets(buf, 1000, ptr_file) != NULL)
		{
			int counter = 0;
			char* ptr = buf;
			char* word = buf;
			int idx[4];
			switch (state)
			{
			case STATE_READ:
				if (strncmp(buf, "volumeelements", 14) == 0)
				{
					state = STATE_VOL;
				}
				else if (strncmp(buf, "points", 6) == 0)
				{
					state = STATE_POINTS;
				}
				break;
			case STATE_VOL:
			{
				int n = atoi(buf);
				state = STATE_ELEM;
			}
			break;
			case STATE_ELEM:
				while (ptr = strchr(word, ' '))
				{
					*ptr = '\0';
					counter++;
					int x = atoi(word);
					if (counter > 2)
					{
						idx[counter - 3] = x - 1;
					}
					word = ptr + 1;
				}
				idx[3] = atoi(word) - 1;
				if (counter == 0)
					state = STATE_READ;
				else
				{
					Tet tet;
					tet.idx[0] = idx[0];
					tet.idx[1] = idx[1];
					tet.idx[2] = idx[2];
					tet.idx[3] = idx[3];
					tets.push_back(tet);
				}
				break;
			case STATE_POINTS:
			{
				int n = atoi(buf);
				//SetNumNodes(n);
				nodes.resize(n);
				state = STATE_POINT;
			}
			break;
			case STATE_POINT:
				Vector3R pos;
				while (ptr = strchr(word, ' '))
				{
					if (ptr == word)
					{
						word++;
						continue;
					}
					*ptr = '\0';
					real x = (real)atof(word);
					pos[counter] = x;
					counter++;
					word = ptr + 1;
				}
				pos[2] = (real)atof(word);
				if (counter == 0)
					state = STATE_READ;
				else
					nodes.at(nod++).pos = scale * pos;
				break;
			}
		}

		fclose(ptr_file);
		return true;
	}

	bool IO::LoadFromTet1File(const char* path, real scale, const Vector3R& offset, std::vector<Node>& nodes, std::vector<Tet>& tets)
	{
		FILE* ptr_file;
		char buf[1000];

		errno_t ret = fopen_s(&ptr_file, path, "r");
		if (ret != 0)
		{
			printf("Couldn't open file %s\n", path);
			return false;
		}
		const int STATE_READ = 0;
		const int STATE_VOL = 1;
		const int STATE_ELEM = 2;
		const int STATE_POINTS = 3;
		const int STATE_POINT = 4;
		int state = STATE_READ;
		int nod = 0;
		while (fgets(buf, 1000, ptr_file) != NULL)
		{
			char* ptr = buf;
			char* word = buf + 2;

			if (buf[0] == 'v')
			{
				// read node position
				Vector3R pos;
				int counter = 0;
				while (ptr = strchr(word, ' '))
				{
					if (ptr == word)
					{
						word++;
						continue;
					}
					*ptr = '\0';
					real x = (real)atof(word);
					pos[counter] = x;
					counter++;
					word = ptr + 1;
				}
				pos[2] = (real)atof(word);
				std::swap(pos[2], pos[1]); // hack for inverting axes
				nodes.push_back(Node());
				nodes.at(nod++).pos = scale * pos + offset;
			}

			if (buf[0] == 't')
			{
				// read tetrahedral element
				int indx[4];
				int counter = 0;
				while (ptr = strchr(word, ' '))
				{
					*ptr = '\0';

					int x = atoi(word);
					indx[counter] = x;
					counter++;
					word = ptr + 1;
				}
				indx[3] = atoi(word);

				Tet tet;
				tet.idx[0] = indx[0];
				tet.idx[1] = indx[1];
				tet.idx[2] = indx[2];
				tet.idx[3] = indx[3];
				tets.push_back(tet);
			}
		}

		fclose(ptr_file);
		return true;
	}

	bool IO::LoadFromXmlFile(const char* path, std::vector<Node>& nodes, std::vector<Tet>& tets,
		std::vector<int>& fixedNodes, std::vector<uint32>& surfTris, FemConfig& femConfig)
	{
		// extract directory name
		std::string fileName(path);
		std::string dir = fileName.substr(0, fileName.find_last_of("/\\") + 1);

		fixedNodes.clear();
		nodes.clear();
		tets.clear();
		surfTris.clear();

		tinyxml2::XMLDocument doc;
		if (doc.LoadFile(path) != tinyxml2::XML_SUCCESS)
			return false;

		tinyxml2::XMLElement* root = doc.FirstChildElement("solidfem");
		if (!root)
			return false;

		tinyxml2::XMLElement* xModel = root->FirstChildElement("model");
		if (xModel)
		{
			// load the model
			const char* modelName = xModel->Attribute("path");
			if (!modelName)
				return false;
			std::string modelPath = dir + modelName;
			real scale = 1;
			const char* xScale = xModel->Attribute("scale");
			if (xScale)
				scale = atof(xScale);

			const char* type = xModel->Attribute("type");
			if (strcmp(type, "vol") == 0)
			{
				if (!LoadFromVolFile(modelPath.c_str(), scale, nodes, tets))
					return false;
			}
			else if (strcmp(type, "tet1") == 0)
			{
				if (!LoadFromTet1File(modelPath.c_str(), scale, Vector3R(), nodes, tets))
					return false;
			}
			else if (strcmp(type, "feb") == 0)
			{
				std::set<uint32> innerSurface;
				return LoadFromFebFile(modelPath.c_str(), nodes, tets, fixedNodes, surfTris, innerSurface, scale, &femConfig);
				// TODO: bcFlag and bcIndices
			}

			int swapCoords = SCT_SWAP_NONE;
			const char* swap = xModel->Attribute("swap");
			if (swap)
			{
				if (strcmp(swap, "xy") == 0)
					swapCoords = SCT_SWAP_X_AND_Y;
				else if (strcmp(swap, "xz") == 0)
					swapCoords = SCT_SWAP_X_AND_Z;
				else if (strcmp(swap, "yz") == 0)
					swapCoords = SCT_SWAP_Y_AND_Z;
			}

			// post-process nodes
			FEM_SYSTEM::Vector3R initialVel(0.f, 0, 0.f);
			for (uint32 i = 0; i < nodes.size(); i++)
			{
				FEM_SYSTEM::Node& node = nodes[i];
				if (swapCoords == SCT_SWAP_X_AND_Y)
					std::swap(node.pos.x, node.pos.y);
				else if (swapCoords == SCT_SWAP_X_AND_Z)
					std::swap(node.pos.x, node.pos.z);
				else if (swapCoords == SCT_SWAP_Y_AND_Z)
					std::swap(node.pos.y, node.pos.z);
			}

			// load the Dirichlet BCs
			tinyxml2::XMLElement* xDirichlet = xModel->FirstChildElement("dirichlet");
			if (xDirichlet)
			{
				const char* text = xDirichlet->GetText();
				std::string str(text);
				while (!str.empty())
				{
					size_t next;
					int idx = std::stoi(str, &next);
					fixedNodes.push_back(idx);
					if (next >= str.size())
						break;
					if (str[next] == ',')
						next++;
					str = str.substr(next);
				}

				Printf("%d pinned nodes\n", fixedNodes.size());
			}
		}

		tinyxml2::XMLElement* xMaterial = root->FirstChildElement("material");
		if (xMaterial)
		{
			const char* name = xMaterial->Attribute("name");
			if (strcmp(name, "neo-hookean") == 0)
			{
				femConfig.mMaterial = (MaterialModelType)MMT_NEO_HOOKEAN;
			}

			const char* str = xMaterial->Attribute("young");
			if (str)
			{
				femConfig.mYoungsModulus = atof(str);
			}

			str = xMaterial->Attribute("poisson");
			if (str)
			{
				femConfig.mPoissonRatio = atof(str);
			}

			str = xMaterial->Attribute("density");
			if (str)
			{
				femConfig.mDensity = atof(str);
			}
		}

		tinyxml2::XMLElement* xSim = root->FirstChildElement("simulation");
		if (xSim)
		{
			const char* type = xSim->Attribute("type");
			if (strcmp(type, "quasi-static") == 0)
			{
				femConfig.mSimType = ST_QUASI_STATIC;
			}
			else if (strcmp(type, "explicit") == 0)
			{
				femConfig.mSimType = ST_EXPLICIT;
			}

			const char* method = xSim->Attribute("method");
			if (strcmp(method, "nonlinear") == 0)
			{
				femConfig.mType = MT_NONLINEAR_ELASTICITY;
			}
			else if (strcmp(method, "mixed") == 0)
			{
				femConfig.mType = MT_INCOMPRESSIBLE_NONLINEAR_ELASTICITY;
			}

			const char* str = xSim->Attribute("steps");
			if (str)
			{
				int n = atoi(str);
				if (femConfig.mSimType < ST_EXPLICIT)
					femConfig.mForceApplicationStep = 1.0 / n;
				else
					femConfig.mNumSubsteps = n;
			}

			str = xSim->Attribute("gravity");
			if (str)
			{
				femConfig.mGravity = atof(str);
			}

			tinyxml2::XMLElement* xPressure = xSim->FirstChildElement("pressure");
			if (xPressure)
			{
				const char* pressure = xPressure->Attribute("p");
				if (pressure)
				{
					femConfig.mAppliedPressure = atof(pressure);
				}

				const char* text = xPressure->GetText();
				std::string str(text);
				while (!str.empty())
				{
					size_t next;
					int idx;
					try
					{
						idx = std::stoi(str, &next);
						surfTris.push_back(idx);
					}
					catch (...)
					{
						break;
					}
					if (next >= str.size())
						break;
					str = str.substr(next);
				}
			}
		}

		tinyxml2::XMLElement* xSolver = root->FirstChildElement("solver");
		if (xSolver)
		{
			const char* iters = xSolver->Attribute("iterations");
			if (iters)
			{
				femConfig.mOuterIterations = atoi(iters);
			}

			const char* tol = xSolver->Attribute("tolerance");
			if (tol)
			{
				femConfig.mAbsNewtonRsidualThreshold = atof(tol);
			}
		}
		
		return true;
	}
} // namespace FEM_SYSTEM
