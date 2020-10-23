from TorchSimulator import TorchSimulator, TorchSimulatorFast
from utils import read_tetfile
import time

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
#ts_slow_begin = time.time()
#sim2 = TorchSimulator(verts, indices, fixed2, config)
#sim2.solve_grad_desc(1e-5, 100)
#ts_slow_end = time.time()
#print(sim2.nodes[-1])
#print("Time (slow):", ts_slow_end - ts_slow_begin)

ts_fast_begin = time.time()
sim2 = TorchSimulatorFast(verts, indices, fixed2, config)
sim2.solve_grad_desc(1e-5, 100)
ts_fast_end = time.time()
print(sim2.nodes[-1])
print("Time (fast):", ts_fast_end - ts_fast_begin)
