import pysolidfem as psf
import numpy as np
import matplotlib.pyplot as plt
import torch
import scipy.optimize as opt

from mpl_toolkits.mplot3d import Axes3D
from Simulator import Simulator
from TorchSimulator import TorchSimulator
from utils import read_tetfile, read_objfile

from scipy.optimize import minimize

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
    tets = np.array([0, 1, 2, 3], dtype=np.int)
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
        "young": 10000,
        "poisson": 0.49,
        "simtype": 0,
        "substeps": 1,
        'maxiter': 100,
    }

    # tetrahedron
    hw = 0.05
    tets = np.array([0, 1, 2, 3], dtype=np.int)
    nodes = np.array([[0, -hw, hw], [0, -hw, -hw], [0, hw, hw], [2 * hw, -hw, hw]], dtype=np.float)
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
    #sim.solve_grad_desc(2e-3, 100, 1e-12)
    #sim.solve_steep_desc(num_iters=10)
    #sim.solve_conj_grad(2)
    sim.solve_scipy()
    print(sim.nodes)

    # parameter estimation
    sim1 = Simulator(nodes, tets, fixed_nodes, config)
    sim1.mu = 100
    sim1.la = 1000
    sim1.estimate_params(sim.nodes)
    print(sim.mu, sim.la)
    print(sim1.nodes)

def test_cantilever():
    config = {
        "young": 66000,
        "poisson": 0.45,
        "simtype": 2,
        "substeps": 10,
        'maxiter': 100,
    }

    # load the box
    verts, indices = read_tetfile('../Models/box2.tet')
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

    # Python torch simulator
    sim2 = TorchSimulator(verts, indices, fixed2, config, device='cpu')
    for iter in range(0, num_steps):
        sim2.step()
    print(sim2.nodes[-1])
    plot_nodes(sim2.nodes)

def test_cantilever_static():
    config = {
        "young": 66000,
        "poisson": 0.45,
        'density': 1070,
        "simtype": 0,
        "substeps": 10,
        'maxiter': 100,
        'tol': 0.01,
        'mixed': False,
    }

    # load the box
    verts, indices = read_tetfile('../Models/box2.tet')
    fixed = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    #all = range(0, verts.shape[0])
    #free = np.setdiff1d(all, fixed)
    #print('free: ', free)

    num_steps = 10

    # Python simulator
    sim = Simulator(verts, indices, fixed, config)
    #sim.solve_scipy()
    #sim.solve_conj_grad(400)
    #print(sim.nodes[-1])

    # C++ simulator
    simC = psf.NonlinearFEM(indices, verts, fixed, config)
    #simC = psf.NonlinearFEM('cantilever.xml')
    nodes2 = simC.get_nodes()
    simC.step()
    nodes2 = simC.get_nodes()
    simC.save_to_vtk('cantilever.vtk')
    simC.save_to_obj('cantilever.obj')
    print(nodes2[-1])
    plot_nodes(nodes2)

    V, F, boundary = simC.get_boundary_mesh()
    target_set = np.intersect1d(sim.free, boundary)
    print(target_set)
    target_setf = target_set - len(fixed)

    # parameter estimation
    target = nodes2[target_set]

    def inverse_objective(x):
        mu = x[0]
        la = x[1]
        config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
        config["poisson"] = 0.5 * la / (la + mu)
        simC = psf.NonlinearFEM(indices, verts, fixed, config)
        simC.step()
        nodesC = simC.get_nodes()
        delta = (nodesC[target_set] - target).flatten()
        error = 0.5 * np.dot(delta, delta)
        return error

    # TODO: avoid another simulator pass
    def jacobian(x):
        mu = x[0]
        la = x[1]
        config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
        config["poisson"] = 0.5 * la / (la + mu)
        simC = psf.NonlinearFEM(indices, verts, fixed, config)
        simC.step() # simulate with current positions and params
        nodesC = simC.get_nodes()
        delta = (nodesC[target_set] - target).flatten()
        H = simC.get_hessian()
        simC.compute_force_param_grads()
        f_mu = simC.get_force_mu_grad()
        f_lambda = simC.get_force_lambda_grad()
        grad = np.linalg.solve(H, np.column_stack((f_mu, f_lambda)))
        grad_mu = grad[:,0].reshape(-1, 3)[target_setf].flatten()
        grad_la = grad[:,1].reshape(-1, 3)[target_setf].flatten()
        grad_restricted = np.column_stack((grad_mu, grad_la))
        J = np.matmul(np.transpose(grad_restricted), delta)
        return J

    def residual(x):
        mu = x[0]
        la = x[1]
        config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
        config["poisson"] = 0.5 * la / (la + mu)
        simC = psf.NonlinearFEM(indices, verts, fixed, config)
        simC.step() # simulate with current positions and params
        nodesC = simC.get_nodes()
        delta = (nodesC[target_set] - target).flatten()
        return delta

    def res_jac(x):
        mu = x[0]
        la = x[1]
        config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
        config["poisson"] = 0.5 * la / (la + mu)
        simC = psf.NonlinearFEM(indices, verts, fixed, config)
        simC.step() # simulate with current positions and params
        nodesC = simC.get_nodes()
        delta = (nodesC[target_set] - target).flatten()
        H = simC.get_hessian()
        simC.compute_force_param_grads()
        f_mu = simC.get_force_mu_grad()
        f_lambda = simC.get_force_lambda_grad()
        grad = np.linalg.solve(H, np.column_stack((f_mu, f_lambda)))
        grad_mu = grad[:,0].reshape(-1, 3)[target_setf].flatten()
        grad_la = grad[:,1].reshape(-1, 3)[target_setf].flatten()
        grad_restricted = np.column_stack((grad_mu, grad_la))
        return grad_restricted

    # plot the loss function w.r.t. mu
    #mu_min = 10000
    #mu_max = 30000
    #steps = 1000
    #loss = np.empty(steps, dtype=np.double)
    #mus = np.empty(steps, dtype=np.double)
    #for i in range(0, steps):
    #    mus[i] = mu_min + (mu_max - mu_min) * (i + 1) / steps
    #    x = [mus[i], sim.la]
    #    loss[i] = inverse_objective(x)
    #fig, ax = plt.subplots()
    #ax.plot(mus, loss)
    #fig.savefig("loss.png")

    # plot the loss function w.r.t. lambda
    #la_min = 100000
    #la_max = 300000
    #steps = 5
    #loss = np.empty(steps, dtype=np.double)
    #las = np.empty(steps, dtype=np.double)
    #for i in range(0, steps):
    #    las[i] = la_min + (la_max - la_min) * (i + 1) / steps
    #    x = [sim.mu, las[i]]
    #    loss[i] = inverse_objective(x)
    #fig, ax = plt.subplots()
    #ax.plot(las, loss)
    #fig.savefig("loss_lambda.png")

    mu = 20000;
    la = 200000;
    x0 = [mu, la]
    #sol = minimize(inverse_objective, x0, method='Nelder-Mead')
    #sol = minimize(inverse_objective, x0, method='Powell')
    #sol = minimize(inverse_objective, x0, method="BFGS", options={'gtol': 1e-12}, jac=jacobian)
    sol = opt.least_squares(residual, x0, method='lm', gtol=1e-12, jac=res_jac)
    #sol = opt.dual_annealing(inverse_objective, bounds=((10000, 30000), (200000, 210000)))
    print(sol)
    mu = sol.x[0]
    la = sol.x[1]
    print(mu, la)
    print(sim.mu, sim.la)

    # Python torch simulator
    #sim2 = TorchSimulator(verts, indices, fixed2, config)
    #sim2.solve_grad_desc(3e-5, 70000, 1e-15, 0.1)
    #print(sim2.nodes[-1])

def test_cantilever_cable():
    path = '../Models/cantilever_cable.xml'
    # C++ simulator
    simC = psf.NonlinearFEM(path)
    simC.set_cable_actuation(0, 0.5)
    simC.step()
    nodes2 = simC.get_nodes()
    simC.save_to_vtk('cantilever_cable.vtk')
    simC.save_to_obj('cantilever_cable.obj')

    # Python simulator
    N, T, fixed = psf.load_from_xml(path)
    #config = {
    #    "young": 66000,
    #    "poisson": 0.45,
    #    'density': 1070,
    #    "simtype": 0,
    #    "substeps": 10,
    #}
    #sim = Simulator(N, T.flatten(), fixed, config)

    ## compute target set
    #V, F, boundary = simC.get_boundary_mesh()
    #target_set = np.intersect1d(sim.free, boundary)
    #target_setf = target_set - len(fixed)

    target_set = 76

    # cable control
    target = nodes2[target_set]
    target[0] += 0.1
    print(target)

    def inverse_objective(x):
        alpha = x[0]
        simC = psf.NonlinearFEM(path)
        simC.set_cable_actuation(0, alpha)
        simC.step()
        nodesC = simC.get_nodes()
        delta = (nodesC[target_set] - target).flatten()
        error = 0.5 * np.dot(delta, delta)
        return error

    alpha = 1;
    x0 = [alpha]
    sol = minimize(inverse_objective, x0, method='Nelder-Mead')
    print(sol)
    alpha = sol.x[0]
    print(alpha)

    # simulate again
    simC = psf.NonlinearFEM(path)
    simC.set_cable_actuation(0, alpha)
    simC.step()
    simC.save_to_vtk('cantilever_control.vtk')
    simC.save_to_obj('cantilever_control.obj')

    # plot the mesh
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')                
    #pos = nodesC.flatten()
    #x = pos[0::3]
    #y = pos[1::3]
    #z = pos[2::3]
    #ax.plot_trisurf(x, z, y, linewidth=0.2, antialiased=True, triangles=T)
    #plt.savefig('tetmesh.png')


def test_hammerbot():
    config = {
        "young": 90000,
        "poisson": 0.4999,
        "simtype": 0,
        "substeps": 10,
        'maxiter': 100,
        "mixed" : True,
        'tol': 0.01,
    }

    #mu = 17421.6
    #la = 5.38971e10
    #mu = 22934.3
    #la = 227048.6
    #config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
    #config["poisson"] = 0.5 * la / (la + mu)
    #print(config)

    # create the simulator from file
    path = '../Models/hammerbot.xml'
    fem = psf.NonlinearFEM(path, config)
    fem.step()
    nodes = fem.get_nodes()
    fem.save_to_vtk("../Models/hammerbot.vtk")
    fem.save_to_obj('hammerbot.obj')
    
    N, T, F = psf.load_from_xml(path)
    #print(N)
    hammerbot = Simulator(N, T.flatten(), F, config)
    #nodes4 = hammerbot.step()
    #plot_nodes(nodes4)

    # parameter estimation
    #target_set = range(0, nodes.shape[0])
    #print(target_set)

    V, F, boundary = fem.get_visual_mesh()
    #target = V
    V1, F1 = read_objfile("../Models/hammer_target.obj")
    target = 0.01 * V1
    #target = nodes

    def inverse_objective(x):
        mu = x[0]
        la = x[1]
        config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
        config["poisson"] = 0.5 * la / (la + mu)
        simC = psf.NonlinearFEM(path, config)
        simC.step() # simulate with current positions and params
        V, _, _ = simC.get_visual_mesh()
        delta = (V - target).flatten()
        error = np.dot(delta, delta)
        #print('err {:E}'.format(error))
        return error

    #def residual(x):
    #    mu = x[0]
    #    la = x[1]
    #    config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
    #    config["poisson"] = 0.5 * la / (la + mu)
    #    simC = psf.NonlinearFEM(path, config)
    #    simC.step() # simulate with current positions and params
    #    nodesC = simC.get_nodes()
    #    delta = (nodesC - target).flatten()
    #    return delta

    #def res_jac(x):
    #    mu = x[0]
    #    la = x[1]
    #    config["young"] = mu * (3.0 * la + 2.0 * mu) / (mu + la)
    #    config["poisson"] = 0.5 * la / (la + mu)
    #    simC = psf.NonlinearFEM(path, config)
    #    simC.step() # simulate with current positions and params
    #    nodesC = simC.get_nodes()
    #    H = simC.get_hessian()
    #    simC.compute_force_param_grads()
    #    f_mu = simC.get_force_mu_grad()
    #    f_lambda = simC.get_force_lambda_grad()
    #    grad = np.linalg.solve(H, np.column_stack((f_mu, f_lambda)))
    #    return grad

    mu = 20000;
    la = 1e7;
    x0 = [mu, la]
    sol = minimize(inverse_objective, x0, method='Nelder-Mead')
    #sol = opt.least_squares(residual, x0, method='dogbox', gtol=1e-12, xtol=None, ftol=None) #, jac=res_jac)
    print(sol.message)
    mu = sol.x[0]
    la = sol.x[1]
    print(mu, la)      

    # plot the loss function w.r.t. mu
    #mu_min = 10000
    #mu_max = 50000
    #steps = 1000
    #loss = np.empty(steps, dtype=np.double)
    #mus = np.empty(steps, dtype=np.double)
    #for i in range(0, steps):
    #    mus[i] = mu_min + (mu_max - mu_min) * (i + 1) / steps
    #    x = [mus[i], fem.get_lame_lambda()]
    #    loss[i] = inverse_objective(x)
    #fig, ax = plt.subplots()
    #ax.plot(mus, loss)
    #fig.savefig("loss_hb_mu.png")

    # plot the loss function w.r.t. lambda
    #la_min = 100000
    #la_max = 200000
    #steps = 100
    #loss = np.empty(steps, dtype=np.double)
    #las = np.empty(steps, dtype=np.double)
    #for i in range(0, steps):
    #    las[i] = la_min + (la_max - la_min) * (i + 1) / steps
    #    x = [fem.get_shear_modulus(), las[i]]
    #    loss[i] = inverse_objective(x)
    #fig, ax = plt.subplots()
    #ax.plot(las, loss)
    #fig.savefig("loss_hb_lambda.png")

    # load the hammerbot
    #verts, indices = read_tetfile('hammerbot_fine.tet')
    #fixed2 = [261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 1308, 1309, 1310, 1311, 1312, 1313, 1322, 1323, 1324, 1325, 1334, 1335, 1336, 1337, 1357, 1436, 1492, 1496, 1523, 1564, 1600, 1602, 1656, 1663, 1682, 1800, 1804, 1894, 2168, 2240, 2260]
    #sim2 = TorchSimulator(0.01 * verts, indices, fixed2, config, 'cuda')
    #sim2.solve_grad_desc(2e-4, 70000, 1e-15, 0.1)
    #plot_nodes(sim2.nodes)

def test_hammerbot_cable():
    # create the simulator from file
    path = '../Models/hammerbot_cable.xml'
    fem = psf.NonlinearFEM(path)
    fem.set_cable_actuation(0, 0.82)
    fem.simulate_static()
    nodes = fem.get_nodes()
    fem.save_to_vtk("../Models/hammerbot_cable.vtk")
    fem.save_to_obj('../Models/hammerbot_cable.obj')

torch.set_printoptions(precision=8)
#test_tetrahedron_static()
#test_cantilever_static()
#test_cantilever_cable()
test_hammerbot_cable()

# explicit integration
#fem = psf.NonlinearFEM('hammerbot_explicit.xml')
#fem.save_to_vtk("hammerbot0.vtk")
#for i in range(0, 10):
#    fem.step(dt)
#    fem.save_to_vtk("hammerbot{0}.vtk".format(i + 1))
