#include <pybind11/pybind11.h>
//TODO: use a relative path for the include dir
#include <Include/FemPhysicsMatrixFree.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using namespace FEM_SYSTEM;

class PyNonlinearFEM
{
public:
	PyNonlinearFEM(py::array_t<int> tets, py::array_t<double> nodes, py::array_t<int> fixed_nodes)
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

		// create the FEM object
		FemConfig config; // default config
		mPhys.reset(new FemPhysicsMatrixFree(nTets, nNodes, config));
	}

	void Step()
	{
		mPhys->Step(0.016);
	}

	py::array_t<double> GetNodes() const
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

private:
	std::unique_ptr< FemPhysicsMatrixFree> mPhys;
};

PYBIND11_MODULE(pysolidfem, m) {	
	py::class_<PyNonlinearFEM>(m, "NonlinearFEM")
		.def(py::init<py::array_t<int>, py::array_t<double>, py::array_t<int>>())
		.def("step", &PyNonlinearFEM::Step)
		.def("get_nodes", &PyNonlinearFEM::GetNodes);
}
