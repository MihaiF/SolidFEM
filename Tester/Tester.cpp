// Tester.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "../Include/FemPhysicsLinearElasticity.h"

using namespace FEM_SYSTEM;

int main()
{
    std::cout << "Hello World!\n"; 
	FemPhysicsBase* femPhysics = nullptr;
	
	std::vector<Tetrahedron> tets;
	std::vector<Node> nodes;
	FemConfig config;
	//femPhysics = new FemPhysicsLinearElasticity(tets, nodes, config);
	//delete femPhysics;

	tets.resize(1);
	nodes.resize(4);
	//femPhysics = new FemPhysicsLinearElasticity(tets, nodes, config);
	//delete femPhysics;

	tets[0].i[0] = 0;
	tets[0].i[1] = 1;
	tets[0].i[2] = 2;
	tets[0].i[3] = 3;
	femPhysics = new FemPhysicsLinearElasticity(tets, nodes, config);
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
