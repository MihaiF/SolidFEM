import numpy as np
import torch
from quat import Quaternion as Q

def read_tetfile(tetfile):
    verts = np.empty((0, 3), dtype=np.float)
    indices = np.empty((0, 4), dtype=np.int64)
    with open(tetfile, 'r') as fh:
        lines = fh.readlines()
        for l in lines:
            tokens = l.split()
            if len(tokens) == 0:
                continue
            if tokens[0] == 'v':
                verts = np.append(verts, [list(map(lambda x: float(x.strip()), tokens[1:4]))], axis=0)
            if tokens[0] == 't':
                indices = np.append(indices, list(map(lambda x: int(x.strip()), tokens[1:5])))
    return verts, indices

def read_objfile(objfile):
    verts = np.empty((0, 3), dtype=np.float)
    indices = np.empty((0, 3), dtype=np.int64)
    with open(objfile, 'r') as fh:
        lines = fh.readlines()
        for l in lines:
            tokens = l.split()
            if len(tokens) == 0:
                continue
            if tokens[0] == 'v':
                verts = np.append(verts, [list(map(lambda x: float(x.strip()), tokens[1:4]))], axis=0)
            if tokens[0] == 'f':
                indices = np.append(indices, [list(map(lambda x: int(x.strip().split('/')[0]) - 1, tokens[1:4]))], axis=0)

    return verts, indices

def get_index_map(tetverts, objverts):
    correspondance = {}
    for i, vo in enumerate(objverts):
        for j, vt in enumerate(tetverts):
            if np.linalg.norm(vo - vt) == 0:
                correspondance[i] = j
                break
    return correspondance

def construct_visual_indices(index_map, visual_indices):
    corrected_visual_indices = np.empty((0, 3), dtype=np.int64)
    for vi in visual_indices:
        corrected_visual_indices = np.append(corrected_visual_indices, [[index_map[vi[0]],
                                                                         index_map[vi[1]],
                                                                         index_map[vi[2]]]], axis=0)
    return corrected_visual_indices


def create_vertex_mask(bounding_boxes, vertices):
    mask = np.empty((0, 1), dtype=np.int64)
    for bb in bounding_boxes:
        for vi, v in enumerate(vertices):
            if v[0] >= bb[0, 0] and v[0] <= bb[0, 1] and\
               v[1] >= bb[1, 0] and v[1] <= bb[1, 1] and\
               v[2] >= bb[2, 0] and v[2] <= bb[2, 1]:
                mask = np.append(mask, vi)
    return mask


def interpolate_poses(pose_0, pose_1, starting_position, steps, dt, duration, mask, device='cuda'):
    """
    Creates a set of interpolated frames, based on a set of starting positions and
    two "poses", given as (position, rotation) tuples.
    Returns the positional dt between the two.

    @param pose_0: (position, rotation) pair for the starting position.
    @param pose_1: (position, rotation) pair for the end of the motion.
    @param starting_position: The starting position of the nodes.
    @param steps: The number of intermediate steps between pose_0 and pose_1.
    @param dt: The time between each step.
    @param duration: The duration of the movement.
    @param mask: Mask of the vertices that are affected by pose changes.
    @param device: Torch adapter.

    @return: Positional updates needed to go from pose_0 to pose_1 in a given number of steps.
    """
    key_frames = torch.zeros((steps, mask.shape[0], 3), dtype=torch.float, device=device)
    key_frame_updates = torch.zeros((steps, mask.shape[0], 3), dtype=torch.float, device=device)
    if mask.shape[0] == 0:
        return key_frame_updates
    for i in range(0, steps):
        q_slerped = Q.SLERP(pose_0[1].q, pose_1[1].q, t=i*dt / duration)
        p_slerped = (pose_0[0] + pose_1[0]) * (i*dt / duration)
        next_frame = torch.mm(starting_position[mask] - torch.mean(starting_position[mask], dim=0),
                              torch.from_numpy(q_slerped.get_rotation_matrix()).float().to(device)) +\
                      torch.mean(starting_position[mask], dim=0)
        next_frame = next_frame + p_slerped
        if i != 0:
            key_frame_updates[i] = next_frame - key_frames[i-1]
        key_frames[i, :, :] = next_frame
    return key_frame_updates


