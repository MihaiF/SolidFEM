import pysolidfem as psf
import numpy as np

print('Hello!')

hw = -0.05
tets = np.array([[0, 1, 2, 3]], dtype=np.int)
nodes = np.array([[0, -hw, hw], [0, -hw, -hw], [0, hw, hw], [2 * hw, -hw, hw]], dtype=np.double)
fem = psf.NonlinearFEM(tets, nodes)
fem.step()

