#include "PyNonlinearFEM.h"

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

py::tuple PyNonlinearFEM::GetBoundaryMesh() const
{
	const Mesh& mesh = mBody.GetBoundaryMesh();
	int numVerts = (int)mesh.vertices.size();
	py::array_t<double> pyVerts = py::array_t<double>(numVerts * 3);
	{
		// allocate the buffer
		py::buffer_info buf = pyVerts.request();
		double* ptr = (double*)buf.ptr;
		for (int i = 0; i < numVerts; i++)
		{
			for (int j = 0; j < 3; j++)
				ptr[i * 3 + j] = mesh.vertices[i][j];
		}
		pyVerts.resize({ numVerts, 3 });
	}

	py::list pyFaces = py::cast(mesh.indices);
	py::list pyNodes = py::cast(mBody.GetBoundaryNodes());
	return py::make_tuple(pyVerts, pyFaces, pyNodes);
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
