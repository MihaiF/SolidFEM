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
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "PyNonlinearFEM.h"

namespace py = pybind11;

py::tuple PyLoadFromXml(py::str path)
{
	std::cout << "Running in " << GetCurrentWorkingDir() << std::endl;

	// load the setup file
	std::vector<Node> nodes;
	std::vector<Tet> tets;
	std::vector<int> fixedNodes;
	std::vector<uint32> surfTris;
	FemConfig config; // default config
	std::string visualPath;
	float scale;
	std::vector<CableDescriptor> cables;
	IO::LoadFromXmlFile(std::string(path).c_str(), nodes, tets, fixedNodes, surfTris, config, scale, visualPath, cables);

	int numNodes = (int)nodes.size();
	py::array_t<double> pyNodes = py::array_t<double>(numNodes * 3);
	{
		// allocate the buffer
		py::buffer_info buf = pyNodes.request();
		double* ptr = (double*)buf.ptr;
		for (int i = 0; i < numNodes; i++)
		{
			for (int j = 0; j < 3; j++)
				ptr[i * 3 + j] = nodes[i].pos[j];
		}
		pyNodes.resize({ numNodes, 3 });
	}

	int numTets = (int)tets.size();
	py::array_t<int> pyTets = py::array_t<int>(numTets * 4);
	{
		py::buffer_info buf = pyTets.request();
		int* ptr = (int*)buf.ptr;
		for (int i = 0; i < numTets; i++)
		{
			for (int j = 0; j < 4; j++)
				ptr[i * 4 + j] = tets[i].idx[j];
		}
		pyTets.resize({ numTets, 4 });
	}

	py::list pyFixed = py::cast(fixedNodes);
	return py::make_tuple(pyNodes, pyTets, pyFixed);
}


PYBIND11_MODULE(pysolidfem, m) {
	m.def("load_from_xml", &PyLoadFromXml);
	py::class_<PyNonlinearFEM>(m, "NonlinearFEM")
		.def(py::init<py::array_t<int>, py::array_t<double>, py::array_t<int>, py::dict>())
		.def(py::init<py::str>())
		.def(py::init<py::str, py::dict>())
		.def("step", &PyNonlinearFEM::Step, py::arg("dt") = DT)
		.def("simulate_static", &PyNonlinearFEM::SimulateStatic)
		.def("get_nodes", &PyNonlinearFEM::GetNodes)
		.def("save_to_vtk", &PyNonlinearFEM::SaveToVTK)
		.def("save_to_obj", &PyNonlinearFEM::SaveToOBJ)
		.def("get_boundary_mesh", &PyNonlinearFEM::GetBoundaryMesh)
		.def("get_visual_mesh", &PyNonlinearFEM::GetVisualMesh)
		.def("set_lame_params", &PyNonlinearFEM::SetLameParams)
		.def("get_shear_modulus", &PyNonlinearFEM::GetShearModulus)
		.def("get_lame_lambda", &PyNonlinearFEM::GetLameLambda)
		.def("get_hessian", &PyNonlinearFEM::GetHessian)
		.def("compute_force_param_grads", &PyNonlinearFEM::ComputeForceParamGrads)
		.def("get_force_mu_grad", &PyNonlinearFEM::GetForceMuGrad)
		.def("get_force_lambda_grad", &PyNonlinearFEM::GetForceLambdaGrad)
		.def("get_force_rho_grad", &PyNonlinearFEM::GetForceRhoGrad)
		.def("set_cable_actuation", &PyNonlinearFEM::SetCableActuation);
}
