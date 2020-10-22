import pysolidfem as psf
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from Simulator import Simulator
from utils import read_tetfile

def plot_nodes(nodes):
    # plot the mesh   
    fig = plt.figure()             
    ax = fig.add_subplot(111, projection='3d')
    pos = nodes.flatten()
    x = pos[0::3]
    y = pos[1::3]
    z = pos[2::3]
    ax.scatter(x, z, y)
    plt.savefig('plot.png')

# print something
print('Hello!')

# simulation types
# 0 STATIC,
# 1 QUASI_STATIC,
# 2 EXPLICIT,
# 3 IMPLICIT,

def test_tetrahedron():
    config = {
        "young": 66000,
        "poisson": 0.45,
        "simtype": 2,
        "substeps": 1,
        'maxiter': 100,
    }

    # tetrahedron
    hw = 0.05
    tets = np.array([[0, 1, 2, 3]], dtype=np.int)
    nodes = np.array([[0, -hw, hw], [0, -hw, -hw], [0, hw, hw], [2 * hw, -hw, hw]], dtype=np.double)
    fixed_nodes = np.array([0, 1, 2])

    # create the simulator
    fem = psf.NonlinearFEM(tets, nodes, fixed_nodes, config)

    dt = 0.016
    numSteps = 2
    fem.save_to_vtk('tetrahedron0.vtk')
    for i in range(0, numSteps):
        fem.step(dt)
        nodes1 = fem.get_nodes()
        print(nodes1)
        fem.save_to_vtk('tetrahedron{0}.vtk'.format(i+1))

    # Python simulator
    sim = Simulator(nodes, tets, fixed_nodes, config)
    nodes3 = sim.step_explicit()
    nodes3 = sim.step_explicit()
    print(nodes3)

def test_tetrahedron_static():
    config = {
        "young": 66000,
        "poisson": 0.45,
        "simtype": 0,
        "substeps": 1,
        'maxiter': 100,
    }

    # tetrahedron
    hw = 0.05
    tets = np.array([[0, 1, 2, 3]], dtype=np.int)
    nodes = np.array([[0, -hw, hw], [0, -hw, -hw], [0, hw, hw], [2 * hw, -hw, hw]], dtype=np.double)
    fixed_nodes = np.array([0, 1, 2])

    # create the simulator
    fem = psf.NonlinearFEM(tets, nodes, fixed_nodes, config)

    dt = 0.016
    numSteps = 2
    fem.step(dt)
    nodes1 = fem.get_nodes()
    print(nodes1)
    fem.save_to_vtk('tetrahedron.vtk')

    # Python simulator
    sim = Simulator(nodes, tets, fixed_nodes, config)
    #nodes3 = sim.solve_grad_desc(2e-3, 100, 1e-12)
    nodes3 = sim.solve_scipy()
    print(nodes3)

def test_cantilever():
    config = {
        "young": 66000,
        "poisson": 0.45,
        "simtype": 2,
        "substeps": 10,
        'maxiter': 100,
    }

    # load the box
    verts, indices = read_tetfile('box2.tet')
    fixed2 = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    num_steps = 10

    # C++ simulator
    sim2 = psf.NonlinearFEM(indices, verts, fixed2, config)
    nodes2 = sim2.get_nodes()
    sim2.save_to_vtk('cantilever0.vtk')    
    for iter in range(0, num_steps):
        sim2.step()
        nodes2 = sim2.get_nodes()
        sim2.save_to_vtk('cantilever{0}.vtk'.format(iter + 1))
    print(nodes2[-1])
    plot_nodes(nodes2)

    # Python simulator
    sim2 = Simulator(verts, indices, fixed2, config)
    for iter in range(0, num_steps):
        nodes2 = sim2.step()
    print(nodes2[-1])
    plot_nodes(nodes2)

def test_cantilever_static():
    config = {
        "young": 66000,
        "poisson": 0.45,
        "simtype": 0,
        "substeps": 10,
        'maxiter': 100,
    }

    # load the box
    verts, indices = read_tetfile('box2.tet')
    fixed2 = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    num_steps = 10

    # C++ simulator
    sim2 = psf.NonlinearFEM(indices, verts, fixed2, config)
    nodes2 = sim2.get_nodes()
    sim2.step()
    nodes2 = sim2.get_nodes()
    sim2.save_to_vtk('cantilever.vtk')
    print(nodes2[-1])
    plot_nodes(nodes2)

    # Python simulator
    sim2 = Simulator(verts, indices, fixed2, config)
    nodes2 = sim2.step_static()
    print(nodes2[-1])
    plot_nodes(nodes2)


def test_hammerbot():
    config = {
        "young": 66000,
        "poisson": 0.45,
        "simtype": 0,
        "substeps": 10,
        'maxiter': 100,
    }

    # create the simulator from file
    fem = psf.NonlinearFEM('hammerbot.xml')
    fem.step()
    nodes = fem.get_nodes()
    fem.save_to_vtk("hammerbot.vtk")
    
    #N, T, F = psf.load_from_xml('hammerbot.xml')
    #print(N)
    #hammerbot = Simulator(N, T.flatten(), F, config)
    #nodes4 = hammerbot.step()
    #plot_nodes(nodes4)

test_cantilever_static()

# explicit integration
#fem = psf.NonlinearFEM('hammerbot_explicit.xml')
#fem.save_to_vtk("hammerbot0.vtk")
#for i in range(0, 10):
#    fem.step(dt)
#    fem.save_to_vtk("hammerbot{0}.vtk".format(i + 1))
