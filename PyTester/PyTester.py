import pysolidfem as psf
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

# print something
print('Hello!')

# simulation types
# 0 STATIC,
# 1 QUASI_STATIC,
# 2 EXPLICIT,
# 3 IMPLICIT,


# prepaper the FEM tetrahedron
config = {
    "young": 66000,
    "poisson": 0.4,
    "simtype": 2
}

hw = 0.05
tets = np.array([[0, 1, 2, 3]], dtype=np.int)
nodes = np.array([[0, -hw, hw], [0, -hw, -hw], [0, hw, hw], [2 * hw, -hw, hw]], dtype=np.double)
fixed_nodes = np.array([0, 1, 2]);

# create the simulator
fem = psf.NonlinearFEM(tets, nodes, fixed_nodes, config)

dt = 0.016
fem.save_to_vtk('tetrahedron0.vtk')
for i in range(0, 10):
    fem.step(dt)
    nodes = fem.get_nodes()
    print(nodes)
    fem.save_to_vtk('tetrahedron{0}.vtk'.format(i+1))

# plot the mesh   
fig = plt.figure()             
ax = fig.add_subplot(111, projection='3d')
pos = nodes.flatten()
x = pos[0::3]
y = pos[1::3]
z = pos[2::3]
ax.scatter(x, z, y)
plt.savefig('plot.png')

# create the simulator from file
#fem = psf.NonlinearFEM('hammerbot.xml')
#fem.step()
#nodes = fem.get_nodes()
#fem.save_to_vtk("hammerbot.vtk")

fem = psf.NonlinearFEM('hammerbot_explicit.xml')
fem.save_to_vtk("hammerbot0.vtk")

for i in range(0, 10):
    fem.step(dt)
    fem.save_to_vtk("hammerbot{0}.vtk".format(i + 1))

