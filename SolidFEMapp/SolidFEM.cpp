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

#include <iostream>
#include <Include/FemIO.h>
#include <Include/FemPhysicsMatrixFree.h>
#include <Include/FemPhysicsMixed.h>
#include <Engine/Profiler.h>

FILE* out;

using namespace FEM_SYSTEM;

int main(int argc, char* argv[], char* envp[])
{
	if (argc < 2)
	{
		std::cout << "Please specify a file name" << std::endl;
		return -1;
	}

	// init a default config
	int numSteps = 10;
	bool noLineSearch = false;
	double scale = 1;
	bool optimizer = false;

	FemConfig femConfig;
	femConfig.mYoungsModulus = 100000;
	femConfig.mPoissonRatio = 0.499;
	femConfig.mGravity = 0;
	femConfig.mType = MT_NONLINEAR_ELASTICITY;
	femConfig.mOrder = 1;
	femConfig.mSimType = ST_QUASI_STATIC;
	femConfig.mNumSubsteps = 1;
	femConfig.mCustomConfig = nullptr;
	femConfig.mOuterIterations = 100;
	femConfig.mInnerIterations = 1;
	femConfig.mForceApplicationStep = 1.0 / numSteps;
	femConfig.mHasCollisions = false;
	femConfig.mAbsNewtonRsidualThreshold = 0.1;
	femConfig.mVerbose = true;
	femConfig.mAppliedPressure = 4000;
	femConfig.mZUpAxis = false;

	for (int i = 2; i < argc; i++)
	{
		const char* flag = argv[i];
		if (flag[0] == '-' && flag[1] == 'S')
		{
			flag += 2;
			scale = atof(flag);
		}
	}

	// load the mesh+BCs+config from a FEB file
	std::vector<Node> nodes;
    std::vector<Tet> tets;
    std::vector<int> fixedNodes;
    std::vector<uint32> surfTris;
    std::set<uint32> innerSurface;
	std::vector<uint32> bcIndices;
	int bcFlag;
	Printf("Loading %s...\n", argv[1]);
	size_t len = strlen(argv[1]);
	bool ret = false;
	if (argv[1][len - 3] == 'f' && argv[1][len - 2] == 'e' && argv[1][len - 1] == 'b')
		ret = IO::LoadFromFebFile(argv[1], nodes, tets, fixedNodes, surfTris, innerSurface, scale, &femConfig, &bcFlag, &bcIndices);
	else
		ret = IO::LoadFromXmlFile(argv[1], nodes, tets, fixedNodes, surfTris, femConfig);
	if (ret)
    {
        Printf("Num nodes: %d\n", nodes.size());
        Printf("Num tets: %d\n", tets.size());
		Printf("Num fixed: %d\n", fixedNodes.size());
    }
    else
    {
        std::cout << "Load failed" << std::endl;
        return -1;
    }

	numSteps = (int)ceil(1.0 / femConfig.mForceApplicationStep);

	// FIXME
	for (size_t i = 0; i < nodes.size(); i++)
	{
		nodes[i].pos0 = nodes[i].pos;
	}

	// override or set properties from the command line
	std::string suffix;
	for (int i = 2; i < argc; i++)
	{
		const char* flag = argv[i];
		if (flag[0] == '-' && flag[1] == 'M')
		{
			flag += 2;
			if (strcmp(flag, "mixed") == 0)
			{
				femConfig.mType = MT_INCOMPRESSIBLE_NONLINEAR_ELASTICITY;
				suffix += "_mixed";
			}
			else if (strcmp(flag, "nonlinear") == 0)
			{
				femConfig.mType = MT_NONLINEAR_ELASTICITY;
				suffix += "_nonlin";
			}
		}
		else if (flag[0] == '-' && flag[1] == 'p')
		{
			flag += 2;
			real p = atof(flag);
			femConfig.mAppliedPressure = p;
			suffix = suffix + "_" + flag;
		}
		else if (flag[0] == '-' && flag[1] == 'P')
		{
			flag += 2;
			femConfig.mPoissonRatio = atof(flag);
		}
		else if (flag[0] == '-' && flag[1] == 'N')
		{
			flag += 2;
			numSteps = atoi(flag);
			femConfig.mForceApplicationStep = 1.0 / numSteps;
		}
		else if (flag[0] == '-' && flag[1] == 'I')
		{
			flag += 2;
			femConfig.mOuterIterations = atoi(flag);
		}
		else if (flag[0] == '-' && flag[1] == 'G')
		{
			flag += 2;
			femConfig.mGravity = atof(flag);
		}
		else if (flag[0] == '-' && flag[1] == 'R')
		{
			flag += 2;
			femConfig.mAbsNewtonRsidualThreshold = atof(flag);
		}
		else if (flag[0] == '-' && flag[1] == 'Z')
		{
			femConfig.mZUpAxis = true;
		}
		else if (flag[0] == '-' && flag[1] == 'O')
		{
			optimizer = true;
			suffix += "_opt";
		}
		else if (flag[0] == '-')
		{
			flag += 1;
			if (strcmp(flag, "nols") == 0)
			{
				noLineSearch = true;
				suffix += "_nols";
			}
		}
	}

	std::string path(argv[1]);
	std::string name = path.substr(0, path.size() - 4) + suffix;
	std::string outPath = name + ".txt";
	fopen_s(&out, outPath.c_str(), "wt");

	for (uint32 i = 0; i < fixedNodes.size(); i++)
    {
        int idx = fixedNodes[i];
        nodes[idx].invMass = 0.f;
    }

	Printf("Constructing simulator...\n");
	FemPhysicsBase* femPhysics = nullptr;
	if (femConfig.mType == MT_NONLINEAR_ELASTICITY)
	{
		femConfig.mMaterial = (MaterialModelType)MMT_NEO_HOOKEAN;

		FemPhysicsMatrixFree::Config nonlinConfig;
		nonlinConfig.mSolver = noLineSearch ? NST_NEWTON : NST_NEWTON_LS;
		nonlinConfig.mOptimizer = optimizer;
		femConfig.mCustomConfig = &nonlinConfig;

		femPhysics = new FemPhysicsMatrixFree(tets, nodes, femConfig);
	}
	else if (femConfig.mType == MT_INCOMPRESSIBLE_NONLINEAR_ELASTICITY)
	{
		femConfig.mMaterial = (MaterialModelType)MMT_DISTORTIONAL_OGDEN;
		//femConfig.mInvBulkModulus = 0;
		//femConfig.mPoissonRatio = 0.45;
		FemPhysicsMixed::Config mixedConfig;
		mixedConfig.mSolver = noLineSearch ? NST_NEWTON : NST_NEWTON_LS;
		mixedConfig.mConstraintScale = 1;
		//mixedConfig.mPressureOrder = 0;
		femConfig.mCustomConfig = &mixedConfig;
		femPhysics = new FemPhysicsMixed(tets, nodes, femConfig);
	}

	// add dynamic BCs
	for (size_t i = 0; i < bcIndices.size(); i++)
	{
		femPhysics->AddDirichletBC(bcIndices[i], bcFlag);
	}

	Printf("FEM Config:\n");
	Printf("-----------\n");
	Printf("\tE:\t%g\n", femConfig.mYoungsModulus);
	Printf("\tnu:\t%g\n", femConfig.mPoissonRatio);
	Printf("\tp:\t%g\n", femConfig.mAppliedPressure);
	Printf("\tg:\t%g\n", femConfig.mGravity);
	Printf("\tmethod:\t%d\n", femConfig.mType);
	Printf("\tmat:\t%d\n", femConfig.mMaterial);
	Printf("\tsolver:\t%s\n", noLineSearch ? "Newton" : "Newton + line search");
	Printf("\titers:\t%d\n", femConfig.mOuterIterations);
	Printf("\tsteps:\t%d\n", numSteps);
	Printf("\tres:\t%g\n", femConfig.mAbsNewtonRsidualThreshold);

	if (!surfTris.empty())
    {
        femPhysics->SetBoundaryConditionsSurface(surfTris, -femConfig.mAppliedPressure); // flip the sign for now
    }

	Printf("Solving...\n");
    for (int i = 0; i < numSteps; i++)
    {
		Printf("Step %d\n", i + 1);
        femPhysics->Step(0);
    }

	Printf("Saving VTK file\n");
	std::fstream vtkStream;
	std::string vtkPath = name + ".vtk";
	vtkStream.open(vtkPath, std::fstream::out);
	IO::ExportToVTKHeatMap(vtkStream, femPhysics);
	vtkStream.close();
	
	return 0;
}