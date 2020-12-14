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

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <Include/FemBody.h>
#include <Include/FemPhysicsMatrixFree.h>
#include <Include/FemPhysicsMixed.h>
#include <Include/FemIO.h>

namespace py = pybind11;

using namespace FEM_SYSTEM;

const real DT = 0.016;

class PyNonlinearFEM
{
public:
	PyNonlinearFEM(py::array_t<int> tets, py::array_t<double> nodes, py::array_t<int> fixed_nodes, py::dict config);
	PyNonlinearFEM(py::str path);
	PyNonlinearFEM(py::str path, py::dict config);
	void Step(real dt = DT);
	py::array_t<double> GetNodes() const;
	py::array_t<int> GetTets() const;
	void SaveToVTK(py::str path);
	void SaveToOBJ(py::str path);
	void SetLameParams(real mu, real lambda) { mPhys->SetLameParams(mu, lambda); }
	real GetShearModulus() const { return mPhys->GetShearModulus(); }
	real GetLameLambda() const { return mPhys->GetLameFirstParam(); }
	EigenMatrix GetHessian() const { return EigenMatrix(mPhys->GetHessian()); }
	void ComputeForceParamGrads() { mPhys->GetForceParamGrads(mForceGradMu, mForceGradLambda, mForceRhoLambda); }
	EigenVector GetForceMuGrad() const { return mForceGradMu; }
	EigenVector GetForceLambdaGrad() const { return mForceGradLambda; }
	EigenVector GetForceRhoGrad() const { return mForceGradLambda; }
	py::tuple GetBoundaryMesh() const;
	py::tuple GetVisualMesh() const;
	void SetCableActuation(uint32 cable, real actuation) { mPhys->SetCableActuation(cable, actuation); }

private:
	FemConfig ParseConfig(py::dict config);

private:
	FemBody mBody;
	std::unique_ptr<FemPhysicsBase> mPhys;
	FemPhysicsMatrixFree::Config mNonlinConfig;
	FemPhysicsMixed::Config mMixedConfig;
	bool mUseMixed = false;
	bool mUseOptimizer = false;
	EigenVector mForceGradMu;
	EigenVector mForceGradLambda;
	EigenVector mForceRhoLambda;
};


