from TorchSimulator import TorchSimulator, TorchSimulatorFast
from utils import read_tetfile
import time
import torch

torch.set_printoptions(precision=10)

config = {
    "young": 66000,
    "poisson": 0.45,
    "simtype": 0,
    "substeps": 10,
    'maxiter': 100,
}

num_steps = 100
device = 'cuda'

# load the box
#verts, indices = read_tetfile('box2.tet')
#fixed2 = [0, 1, 2, 3, 4, 5, 6, 7, 8]

# load the hammerbot
verts, indices = read_tetfile('hammerbot_fine.tet')
verts = verts * 0.01
fixed2 = [261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 1308, 1309, 1310, 1311, 1312, 1313, 1322, 1323, 1324, 1325, 1334, 1335, 1336, 1337, 1357, 1436, 1492, 1496, 1523, 1564, 1600, 1602, 1656, 1663, 1682, 1800, 1804, 1894, 2168, 2240, 2260]

# Python torch simulator
print('Creating simulator...')
ts_slow_begin = time.time()
sim2 = TorchSimulator(verts, indices, fixed2, config, device)
ts_slow_end = time.time()
print("Creation time:", ts_slow_end - ts_slow_begin)
print('Simulating...')
ts_slow_begin = time.time()
# for iter in range(0, num_steps):
#     sim2.step_explicit()
#sim2.solve_grad_desc(3e-5, 70000, 1e-15, 0.1) # tuned for the cantilever
sim2.solve_grad_desc(2e-4, 70000, 1e-15, 0.1)
ts_slow_end = time.time()
print(sim2.nodes[-1])
print("Simulation time:", ts_slow_end - ts_slow_begin)

# ts_fast_begin = time.time()
# sim2 = TorchSimulatorFast(verts, indices, fixed2, config, device)
# sim2.solve_grad_desc(1e-5, 100)
# ts_fast_end = time.time()
# print(sim2.nodes[-1])
# print("Time (fast):", ts_fast_end - ts_fast_begin)
