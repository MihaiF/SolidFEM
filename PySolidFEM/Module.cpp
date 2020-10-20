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

#include <pybind11/pybind11.h>
//TODO: use a relative path for the include dir
#include <Include/FemPhysicsMatrixFree.h>
#include <Include/FemIO.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using namespace FEM_SYSTEM;

const real DT = 0.016;

class PyNonlinearFEM
{
public:
	PyNonlinearFEM(py::array_t<int> tets, py::array_t<double> nodes, py::array_t<int> fixed_nodes, py::dict config);
	PyNonlinearFEM(py::str path);
	void Step(real dt = DT) { mPhys->Step(dt); }
	py::array_t<double> GetNodes() const;
	void SaveToVTK(py::str path);

private:
	std::unique_ptr< FemPhysicsMatrixFree> mPhys;
};

PyNonlinearFEM::PyNonlinearFEM(py::array_t<int> tets, py::array_t<double> nodes, py::array_t<int> fixed_nodes, py::dict config)
{
	// read the tets
	py::buffer_info buf = tets.request();
	int rows = (int)buf.shape[0];
	int cols = (int)buf.shape[1]; // should be 4
	if (cols != 4)
	{
		std::cout << "Tets should have 4 columns" << std::endl;
		return;
	}
	// create the native tets
	std::vector<Tet> nTets(rows);
	int* ptr = (int*)buf.ptr;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			nTets[i].idx[j] = ptr[i * cols + j];
		}
	}
	// read the nodes
	buf = nodes.request();
	rows = (int)buf.shape[0];
	cols = (int)buf.shape[1]; // should be 3
	if (cols != 3)
	{
		std::cout << "Nodes should have 3 columns" << std::endl;
		return;
	}
	// create the native nodes
	std::vector<Node> nNodes(rows);
	double* pNodes = (double*)buf.ptr;
	for (int i = 0; i < rows; i++)
	{
		Vector3R pos;
		for (int j = 0; j < cols; j++)
		{
			pos[j] = pNodes[i * cols + j];
		}
		nNodes[i].pos = pos;
		nNodes[i].pos0 = pos;
	}

	// set the fixed nodes
	auto bufFixed = fixed_nodes.request();
	int* ptrFixed = (int*)bufFixed.ptr;
	for (int i = 0; i < fixed_nodes.size(); i++)
	{
		int idx = ptrFixed[i];
		nNodes[idx].invMass = 0;
	}

	// prepare the FEM config
	FemConfig nConfig; // default config
	if (config.contains("young"))
	{
		nConfig.mYoungsModulus = py::cast<real>(config["young"]);
	}
	if (config.contains("poisson"))
	{
		nConfig.mPoissonRatio = py::cast<real>(config["poisson"]);
	}
	if (config.contains("simtype"))
	{ 
		nConfig.mSimType = (SimulationType)py::cast<int>(config["simtype"]);
	}
	nConfig.mMaterial = (MaterialModelType)MMT_NEO_HOOKEAN;
	FemPhysicsMatrixFree::Config nonlinConfig;
	nonlinConfig.mSolver = NST_NEWTON_LS;
	nonlinConfig.mOptimizer = false;
	nConfig.mCustomConfig = &nonlinConfig;

	// create the FEM object
	mPhys.reset(new FemPhysicsMatrixFree(nTets, nNodes, nConfig));
}

PyNonlinearFEM::PyNonlinearFEM(py::str path)
{
	std::cout << "Running in " << GetCurrentWorkingDir() << std::endl;

	// load the setup file
	std::vector<Node> nodes;
	std::vector<Tet> tets;
	std::vector<int> fixedNodes;
	std::vector<uint32> surfTris;
	FemConfig config; // default config
	IO::LoadFromXmlFile(std::string(path).c_str(), nodes, tets, fixedNodes, surfTris, config);

	// FIXME
	for (size_t i = 0; i < nodes.size(); i++)
	{
		nodes[i].pos0 = nodes[i].pos;
	}
	for (uint32 i = 0; i < fixedNodes.size(); i++)
	{
		int idx = fixedNodes[i];
		nodes[idx].invMass = 0.f;
	}

	// create the FEM object
	config.mMaterial = (MaterialModelType)MMT_NEO_HOOKEAN;
	FemPhysicsMatrixFree::Config nonlinConfig;
	nonlinConfig.mSolver = NST_NEWTON_LS;
	nonlinConfig.mOptimizer = false;
	config.mCustomConfig = &nonlinConfig;
	mPhys.reset(new FemPhysicsMatrixFree(tets, nodes, config));
}

py::array_t<double> PyNonlinearFEM::GetNodes() const
{
	// allocate the buffer
	int numNodes = mPhys->GetNumNodes();
	py::array_t<double> result = py::array_t<double>(numNodes * 3);
	py::buffer_info buf = result.request();
	double* ptr = (double*)buf.ptr;
	for (int i = 0; i < numNodes; i++)
	{
		for (int j = 0; j < 3; j++)
			ptr[i * 3 + j] = mPhys->GetDeformedPosition(i)[j];
	}
	result.resize({ numNodes, 3 });
	return result;
}

void PyNonlinearFEM::SaveToVTK(py::str path)
{
	Printf("Saving VTK file\n");
	std::fstream vtkStream;
	std::string vtkPath(path);
	vtkStream.open(vtkPath, std::fstream::out);
	IO::ExportToVTKHeatMap(vtkStream, mPhys.get());
	vtkStream.close();
}

PYBIND11_MODULE(pysolidfem, m) {
	py::class_<PyNonlinearFEM>(m, "NonlinearFEM")
		.def(py::init<py::array_t<int>, py::array_t<double>, py::array_t<int>, py::dict>())
		.def(py::init<py::str>())
		.def("step", &PyNonlinearFEM::Step, py::arg("dt") = DT)
		.def("get_nodes", &PyNonlinearFEM::GetNodes)
		.def("save_to_vtk", &PyNonlinearFEM::SaveToVTK);
}
