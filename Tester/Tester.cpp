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

// Tester.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "../Include/FemPhysicsLinearElasticity.h"
#include "../Include/FemPhysicsMixed.h"
#include "../Include/FemPhysicsMatrixFree.h"

using namespace FEM_SYSTEM;

int main()
{
    std::cout << "Hello World!\n"; 
	FemPhysicsBase* femPhysics = nullptr;
	
	std::vector<Tet> tets;
	std::vector<Node> nodes;
	FemConfig config;
	//femPhysics = new FemPhysicsLinearElasticity(tets, nodes, config);
	//delete femPhysics;

	tets.resize(1);
	nodes.resize(4);
	//femPhysics = new FemPhysicsLinearElasticity(tets, nodes, config);
	//delete femPhysics;

	tets[0].idx[0] = 0;
	tets[0].idx[1] = 1;
	tets[0].idx[2] = 2;
	tets[0].idx[3] = 3;
	femPhysics = new FemPhysicsLinearElasticity(tets, nodes, config);
	delete femPhysics;

	femPhysics = new FemPhysicsMixed(tets, nodes, config);
	delete femPhysics;

	femPhysics = new FemPhysicsMatrixFree(tets, nodes, config);
	delete femPhysics;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
