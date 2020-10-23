from TorchSimulator import TorchSimulator
from utils import read_tetfile

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

# Python torch simulator
sim2 = TorchSimulator(verts, indices, fixed2, config)
sim2.solve_grad_desc(1e-5, 100)
print(sim2.nodes[-1])

