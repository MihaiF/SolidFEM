import pysolidfem as psf
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

print('Hello!')

config = {
    "young": 66000,
    "poisson": 0.4,
}

hw = -0.05
tets = np.array([[0, 1, 2, 3]], dtype=np.int)
nodes = np.array([[0, -hw, hw], [0, -hw, -hw], [0, hw, hw], [2 * hw, -hw, hw]], dtype=np.double)
fixed_nodes = np.array([0, 1, 2]);

fem = psf.NonlinearFEM(tets, nodes, fixed_nodes)
fem.step()
def_nodes = fem.get_nodes()
print(def_nodes)

# plot the mesh   
fig = plt.figure()             
ax = fig.add_subplot(111, projection='3d')
pos = nodes.flatten()
x = pos[0::3]
y = pos[1::3]
z = pos[2::3]
ax.scatter(x, y, z)
plt.savefig('plot.png')


