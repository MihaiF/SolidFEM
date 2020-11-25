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

PyNonlinearFEM::PyNonlinearFEM(py::array_t<int> tets, py::array_t<double> nodes, py::array_t<int> fixed_nodes, py::dict config)
{
	// read the tets
	py::buffer_info buf = tets.request();
	int rows = (int)buf.shape[0];
	int cols = (int)buf.shape[1]; // should be 4 ideally
	if (cols != 4)
	{
		rows = rows / 4;
	}
	// create the native tets
	std::vector<Tet>& nTets = mBody.GetTets();
	nTets.resize(rows);
	int* ptr = (int*)buf.ptr;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			nTets[i].idx[j] = ptr[i * 4 + j];
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
	std::vector<Node>& nNodes = mBody.GetNodes();
	nNodes.resize(rows);
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

	mBody.Prepare(); // fixed nodes are empty for now
	mBody.BuildBoundaryMesh();

	// create the FEM object
	FemConfig nConfig = ParseConfig(config);
	if (mUseMixed)
		mPhys.reset(new FemPhysicsMixed(nTets, nNodes, nConfig));
	else
		mPhys.reset(new FemPhysicsMatrixFree(nTets, nNodes, nConfig));
}

PyNonlinearFEM::PyNonlinearFEM(py::str path)
{
	//std::cout << "Running in " << GetCurrentWorkingDir() << std::endl;

	// load the setup file
	mBody.LoadFromXml(std::string(path).c_str());

	// create the FEM object
	if (!mUseMixed)
	{
		mBody.GetConfig().mMaterial = (MaterialModelType)MMT_NEO_HOOKEAN;
		mNonlinConfig.mSolver = NST_NEWTON_LS;
		mNonlinConfig.mOptimizer = true;
		mBody.GetConfig().mCustomConfig = &mNonlinConfig;
		mPhys.reset(new FemPhysicsMatrixFree(mBody.GetTets(), mBody.GetNodes(), mBody.GetConfig()));
	}
	else
	{
		mBody.GetConfig().mMaterial = (MaterialModelType)MMT_DISTORTIONAL_OGDEN;
		mMixedConfig.mSolver = NST_NEWTON_LS;
		mBody.GetConfig().mCustomConfig = &mMixedConfig;
		mPhys.reset(new FemPhysicsMixed(mBody.GetTets(), mBody.GetNodes(), mBody.GetConfig()));
	}
}

void PyNonlinearFEM::Step(real dt /*= DT*/)
{
	mPhys->Step(dt);
	mBody.UpdateBoundaryMesh(mPhys.get());
	if (mBody.HasVisualMesh())
		mBody.UpdateVisualMesh();
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

py::array_t<int> PyNonlinearFEM::GetTets() const
{
	int numTets = mPhys->GetNumElements();
	py::array_t<int> result = py::array_t<int>(numTets * 4);
	py::buffer_info buf = result.request();
	int* ptr = (int*)buf.ptr;
	for (int i = 0; i < numTets; i++)
	{
		for (int j = 0; j < 4; j++)
			ptr[i * 4 + j] = mPhys->GetGlobalIndex(i, j);
	}
	result.resize({ numTets, 4 });
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

void PyNonlinearFEM::SaveToOBJ(py::str path)
{
	Printf("Saving OBJ file\n");
	std::string objPath(path);
	if (mBody.HasVisualMesh())
		mBody.SaveVisualMesh(objPath.c_str());
	else
		mBody.SaveBoundaryMesh(objPath.c_str());
}

FEM_SYSTEM::FemConfig PyNonlinearFEM::ParseConfig(py::dict config)
{
	FemConfig nConfig; // default config
	if (config.contains("young"))
	{
		nConfig.mYoungsModulus = py::cast<real>(config["young"]);
	}
	if (config.contains("poisson"))
	{
		nConfig.mPoissonRatio = py::cast<real>(config["poisson"]);
	}
	if (config.contains("density"))
	{
		nConfig.mDensity = py::cast<real>(config["density"]);
	}
	if (config.contains("simtype"))
	{
		nConfig.mSimType = (SimulationType)py::cast<int>(config["simtype"]);
	}
	if (config.contains("substeps"))
	{
		nConfig.mNumSubsteps = py::cast<int>(config["substeps"]);
	}
	if (config.contains("maxiter"))
	{
		nConfig.mOuterIterations = py::cast<int>(config["maxiter"]);
	}
	if (config.contains("tol"))
	{
		nConfig.mAbsNewtonRsidualThreshold = py::cast<real>(config["tol"]);
	}
	if (config.contains("mixed"))
	{
		mUseMixed = py::cast<bool>(config["mixed"]);
	}
	if (!mUseMixed)
	{
		nConfig.mMaterial = (MaterialModelType)MMT_NEO_HOOKEAN;
		mNonlinConfig.mSolver = NST_NEWTON_LS;
		mNonlinConfig.mOptimizer = true;
		nConfig.mCustomConfig = &mNonlinConfig;
	}
	else
	{
		nConfig.mMaterial = (MaterialModelType)MMT_DISTORTIONAL_OGDEN;
		mMixedConfig.mSolver = NST_NEWTON_LS;
		nConfig.mCustomConfig = &mMixedConfig;
	}
	return nConfig;
}

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
	IO::LoadFromXmlFile(std::string(path).c_str(), nodes, tets, fixedNodes, surfTris, config, scale, visualPath);

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
		.def("step", &PyNonlinearFEM::Step, py::arg("dt") = DT)
		.def("get_nodes", &PyNonlinearFEM::GetNodes)
		.def("save_to_vtk", &PyNonlinearFEM::SaveToVTK)
		.def("save_to_obj", &PyNonlinearFEM::SaveToOBJ)
		.def("set_lame_params", &PyNonlinearFEM::SetLameParams)
		.def("get_shear_modulus", &PyNonlinearFEM::GetShearModulus)
		.def("get_lame_lambda", &PyNonlinearFEM::GetLameLambda)
		.def("get_hessian", &PyNonlinearFEM::GetHessian)
		.def("compute_force_param_grads", &PyNonlinearFEM::ComputeForceParamGrads)
		.def("get_force_mu_grad", &PyNonlinearFEM::GetForceMuGrad)
		.def("get_force_lambda_grad", &PyNonlinearFEM::GetForceLambdaGrad)
		.def("get_force_rho_grad", &PyNonlinearFEM::GetForceRhoGrad);
}
