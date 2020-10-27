import torch
import numpy as np
import math

from torch_scatter import scatter_add

class TorchSimulator:

    young = 66000
    poisson = 0.45
    density = 1070    

    def __init__(self,
                nodes,
                tets,
                fixed_nodes,
                config,
                device='cpu'
        ):
        """
        @param nodes: Array of node positions
        @param tets: Array of tetrahedra indices
        @param fixed_nodes: Array of indices indicating fixed nodes
        @param config: Dictionary containing configuration
        @param device: CPU or GPU device for torch
        """
        print("CUDA support: {}".format(torch.cuda.is_available()))

        # init members
        self.substeps = 1
        if "substeps" in config:
            self.substeps = config["substeps"]
        self.dt = 0.016 / self.substeps

        self.device = device
        self.gravity = torch.tensor((0, -9.8, 0), dtype=torch.float, device=self.device)

        # init FEM configuration
        self.nodes = torch.from_numpy(nodes).float().to(self.device)
        numTets = math.ceil(tets.shape[0] / 4)
        numNodes = nodes.shape[0]
        self.tets = tets.reshape((numTets, 4))
        self.masses = torch.zeros(numNodes, dtype=torch.float, device=self.device)
        self.vel = torch.zeros(nodes.shape, dtype=torch.float, device=self.device)
        # build shape matrices (batch)
        self.indices = torch.from_numpy(tets).long().to(self.device)
        self.idx0 = self.indices[0::4]
        self.idx1 = self.indices[1::4]
        self.idx2 = self.indices[2::4]
        self.idx3 = self.indices[3::4]
        col1 = self.nodes[self.idx1] - self.nodes[self.idx0]
        col2 = self.nodes[self.idx2] - self.nodes[self.idx0]
        col3 = self.nodes[self.idx3] - self.nodes[self.idx0]
        shape_mats = torch.empty(numTets, 3, 3, device=self.device)
        shape_mats[:,:,0] = col1
        shape_mats[:,:,1] = col2
        shape_mats[:,:,2] = col3
        # compute volumes and masses
        self.vols = torch.det(shape_mats) / 6
        self.inv_mats = torch.inverse(shape_mats)        
        elem_masses = self.vols * self.density / 4
        # average the node masses
        scatter_add(elem_masses, self.idx0, -1, self.masses)
        scatter_add(elem_masses, self.idx1, -1, self.masses)
        scatter_add(elem_masses, self.idx2, -1, self.masses)
        scatter_add(elem_masses, self.idx3, -1, self.masses)
        # compute inverse masses
        self.inv_masses = torch.zeros_like(self.masses, dtype=torch.float, device=self.device)
        free_idx = np.setdiff1d(range(0, nodes.shape[0]), fixed_nodes)
        self.free = torch.from_numpy(free_idx).long().to(self.device)
        self.inv_masses[self.free] = 1.0 / self.masses[self.free]
        # init params
        self.mu = 0.5 * self.young / (1 + self.poisson)
        self.la = self.young * self.poisson / (1 + self.poisson) / (1 - 2 * self.poisson)

    def step(self):
        for iter in range(0, self.substeps):
            self.step_explicit()

    def step_explicit(self):
        # compute forces
        self.compute_def_grads()
        f = self.compute_forces()
        # Semi-implicit (symplectic) Euler integration
        self.vel[self.free] = self.vel[self.free] + self.dt * (self.gravity + self.inv_masses[self.free] * f[self.free])
        self.nodes[self.free] = self.nodes[self.free] + self.dt * self.vel[self.free]

    # gradient descent static solver
    def solve_grad_desc(self, alpha = 1e-2, numIters = 1000, rel_tol=1e-15, abs_tol = 0.1, c1 = 1e-3):
        # gradient descent
        self.compute_def_grads()
        f = self.compute_gradients()
        sqNorm = torch.matmul(f[self.free].view(-1), f[self.free].view(-1))
        energy = self.compute_energy()
        for iter in range(0, numIters):            
            self.nodes[self.free] = self.nodes[self.free] + alpha * f[self.free]
            # update gradients and energy
            self.compute_def_grads()
            f = self.compute_gradients()
            sqNorm = torch.matmul(f[self.free].view(-1), f[self.free].view(-1))
            old_energy = energy
            energy = self.compute_energy()
            rel_err = (energy - old_energy) / energy
            if energy > old_energy + c1 * sqNorm:
                print('Not converging at iteration {}'.format(iter + 1))
                break
            if abs(rel_err) < rel_tol and sqNorm < abs_tol * abs_tol:
                print('Converged after {0} iterations'.format(iter+1))
                break        
        print('energy: {0}, rel change: {1}, grad norm: {2}'.format(energy, rel_err, torch.sqrt(sqNorm)))

    # returns the elastic forces
    def compute_def_grads(self):
        # build shape matrices (batch)
        idx0 = self.indices[0::4]
        idx1 = self.indices[1::4]
        idx2 = self.indices[2::4]
        idx3 = self.indices[3::4]
        col1 = self.nodes[idx1] - self.nodes[idx0]
        col2 = self.nodes[idx2] - self.nodes[idx0]
        col3 = self.nodes[idx3] - self.nodes[idx0]
        numTets = numTets = math.ceil(self.indices.shape[0] / 4)
        shape_mats = torch.empty(numTets, 3, 3, device=self.device)
        shape_mats[:,:,0] = col1
        shape_mats[:,:,1] = col2
        shape_mats[:,:,2] = col3
        #compute deformation gradients (batch)
        self.F = torch.bmm(shape_mats, self.inv_mats)
        self.Finv = torch.inverse(self.F)        
        self.J = torch.det(self.F)
        self.logJ = torch.log(self.J)

    def compute_forces(self):
        # compute the stress tensors (first Piola-Kirchoff)        
        Fit = torch.transpose(self.Finv, 1, 2)
        P = self.mu * (self.F - Fit) + self.la * self.logJ.view(-1,1,1) * Fit
        # compute force matrices
        Pvol = -self.vols.view(-1,1,1) * P
        force_mats = torch.bmm(Pvol, torch.transpose(self.inv_mats, 1, 2))
        
        #numNodes = self.nodes.shape[0]
        f = torch.zeros(self.nodes.shape, dtype=torch.float, device=self.device)
        f0 = -(force_mats[:,:,0] + force_mats[:,:,1] + force_mats[:,:,2])
        scatter_add(f0, self.idx0, 0, f)
        scatter_add(force_mats[:,:,0], self.idx1, 0, f)
        scatter_add(force_mats[:,:,1], self.idx2, 0, f)
        scatter_add(force_mats[:,:,2], self.idx3, 0, f)
        return f

    # returns the negative gradients
    def compute_gradients(self):
        f = self.compute_forces()
        f[self.free] = f[self.free] + self.masses[self.free].unsqueeze(1) * self.gravity
        return f

    def compute_energy(self):
        energy = self.compute_elastic_energy()
        G = self.masses * torch.matmul(self.nodes, self.gravity)
        return energy - G[self.free].sum(-1) # drop the self.free

    def compute_elastic_energy(self):        
        # compute the per element energies
        C = torch.bmm(torch.transpose(self.F, 1, 2), self.F)
        I1 = torch.diagonal(C, dim1=1, dim2=2).sum(-1) # the trace of C
        E = 0.5 * self.mu * (I1 - 3) - self.mu * self.logJ + 0.5 * self.la * self.logJ * self.logJ
        elemE = self.vols * E
        return elemE.sum(-1)

    def compute_elem_energy(self, e):
        F = self.compute_deformation_gradient(e)
        J = torch.det(F)
        if J < 0:
            return float('inf')
        C = torch.matmul(F.T, F)
        I1 = torch.trace(C)
        logJ = math.log(J)
        E = 0.5 * self.mu * (I1 - 3) - self.mu * logJ + 0.5 * self.la * logJ * logJ
        return self.vols[e] * E

    def compute_stress(self, e):
        F = self.compute_deformation_gradient(e)
        # Neo-Hookean
        Fi = torch.inverse(F)
        Fit = Fi.T
        J = torch.det(F)
        P = self.mu * (F - Fit) + self.la * math.log(J) * Fit
        if e == 0:
            print(P)
        return P

    def compute_deformation_gradient(self, e):
        tet = self.tets[e]
        # get the 4 points
        x0 = self.nodes[tet[0]]
        x1 = self.nodes[tet[1]]
        x2 = self.nodes[tet[2]]
        x3 = self.nodes[tet[3]]
        # build the columns of the matrix
        d1 = x1 - x0
        d2 = x2 - x0
        d3 = x3 - x0
        # compute the deformation gradient
        D = torch.cat([d1.view(-1, *d1.shape), d2.view(-1, *d2.shape), d3.view(-1, *d3.shape)], dim=0).T
        F = torch.matmul(D, self.inv_mats[e])
        return F



class TorchSimulatorFast:

    young = 66000
    poisson = 0.45
    density = 1070

    def __init__(self,
                nodes,
                tets,
                fixed_nodes,
                config,
                device='cpu'
        ):
        """
        @param nodes: Array of node positions
        @param tets: Array of tetrahedra indices
        @param fixed_nodes: Array of indices indicating fixed nodes
        @param config: Dictionary containing configuration
        @param device: CPU or GPU device for torch
        """
        print(torch.cuda.is_available())

        # init members
        self.substeps = 1
        if "substeps" in config:
            self.substeps = config["substeps"]
        self.dt = 0.016 / self.substeps

        self.device = device
        self.gravity = torch.tensor((0, -9.8, 0), dtype=torch.float, device=self.device)

        # init FEM configuration
        self.nodes = torch.from_numpy(nodes).float().to(self.device)
        numTets = math.ceil(tets.shape[0] / 4)
        self.tets = tets.reshape((numTets, 4))
        self.masses = torch.zeros((nodes.shape[0], 1), dtype=torch.float, device=self.device)
        self.vel = torch.zeros(nodes.shape, dtype=torch.float, device=self.device)
        # compute volumes and masses
        self.volumes = []
        self.jacobians = []
        self.torch_jacobians = torch.zeros((self.tets.shape[0], 3, 3))
        self.torch_volumes = torch.zeros(self.tets.shape[0])
        for e, tet in enumerate(self.tets):
            # get the 4 points
            x0 = self.nodes[tet[0]]
            x1 = self.nodes[tet[1]]
            x2 = self.nodes[tet[2]]
            x3 = self.nodes[tet[3]]
            # build the columns of the matrix
            d1 = x1 - x0
            d2 = x2 - x0
            d3 = x3 - x0
            # store the volume and the tranpose inverse
            D = torch.cat([d1.view(-1, *d1.shape), d2.view(-1, *d2.shape), d3.view(-1, *d3.shape)], dim=0).T
            vol = torch.det(D) / 6.0
            self.volumes.append(vol)
            self.torch_volumes[e] = vol
            Dinv = torch.inverse(D)
            self.jacobians.append(Dinv.T)
            self.torch_jacobians[e, :, :] = Dinv.T
            # compute the element mass
            mass = 0.25 * vol * self.density
            self.masses[tet[0]] += mass;
            self.masses[tet[1]] += mass;
            self.masses[tet[2]] += mass;
            self.masses[tet[3]] += mass;
        # compute inverse masses
        self.inv_masses = torch.zeros_like(self.masses, dtype=torch.float, device=self.device)
        self.free = np.setdiff1d(range(0, nodes.shape[0]), fixed_nodes)
        self.inv_masses[self.free] = 1.0 / self.masses[self.free]
        # init params
        self.mu = 0.5 * self.young / (1 + self.poisson)
        self.la = self.young * self.poisson / (1 + self.poisson) / (1 - 2 * self.poisson)

    def step(self):
        for iter in range(0, self.substeps):
            self.step_explicit()

    def step_explicit(self):
        # compute forces
        f = self.compute_forces()
        # Semi-implicit (symplectic) Euler integration
        self.vel[self.free] = self.vel[self.free] + self.dt * (self.gravity + self.inv_masses[self.free] * f[self.free])
        self.nodes[self.free] = self.nodes[self.free] + self.dt * self.vel[self.free]

    # gradient descent static solver
    def solve_grad_desc(self, alpha = 1e-2, numIters = 1000, rel_tol=1e-9):
        # gradient descent
        energy = self.compute_energy()
        for iter in range(0, numIters):
            f = self.compute_gradients()
            self.nodes[self.free] = self.nodes[self.free] + alpha * f[self.free]
            #print(energy)
            old_energy = energy;
            energy = self.compute_energy()
            rel_err = abs((energy - old_energy) / energy)
            #print(rel_err)
            if rel_err < rel_tol:
                print('Converged after {0} iterations'.format(iter+1))
                break
            if energy > old_energy:
                print('Not converging')
                break
        print(energy)

    def compute_forces_batched(self):
        P = self.compute_stress_batched()
        H_ = -self.torch_volumes[:, None, None] * torch.bmm(P_, self.torch_jacobians)
        # TODO: implement f


    # returns the elastic forces
    def compute_forces(self):
        f = torch.zeros(self.nodes.shape, dtype=torch.float, device=self.device)
        for e, tet in enumerate(self.tets):
            P = self.compute_stress(e)
            H = -self.volumes[e] * torch.matmul(P, self.jacobians[e])
            for j in range(1,4):
                force = H[:, j - 1]
                f[tet[j]] += force
                f[tet[0]] -= force
        return f

    # returns the negative gradients
    def compute_gradients(self):
        f = self.compute_forces()
        # Inverse masses mask that will filter away all elements with an inverse mass of 0
        # TODO: Make this mask a part of the simulator as it doesn't change between
        # iterations.
        inv_masses = torch.where(self.inv_masses == 0,
                                 torch.tensor(0.),
                                 torch.tensor(1.)).to(self.device)
        f = self.masses * self.gravity * inv_masses
        return f

    def compute_energy(self):
        energy = self.compute_elastic_energy()
        # Inverse masses mask that will filter away all elements with an inverse mass of 0
        # TODO: Make this mask a part of the simulator as it doesn't change between
        # iterations.
        inv_masses = torch.where(self.inv_masses == 0,
                                 torch.tensor(0.),
                                 torch.tensor(1.)).to(self.device).T
        energy = -torch.matmul(self.nodes[:, 1],
                               (self.gravity[1] * self.masses[:, 0] * inv_masses).T)
        return energy

    def compute_elastic_energy_batched(self):
        F = self.compute_deformation_gradient_batched()
        J = torch.det(F)
        if any(J < 0):
            return float('inf')
        C = torch.bmm(torch.transpose(F, -1, -2), F)
        # TODO: Optimize this part - I couldn't find a map/apply function in torch...
        I1 = torch.zeros(C.shape[0], dtype=torch.float, device=self.device)
        for i in range(self.tets.shape[0]):
            I1[i] = torch.trace(C[i])
        logJ = torch.log(J)
        E = 0.5 * self.mu * (I1 - 3) - self.mu * logJ + 0.5 * self.la * logJ * logJ

        energy = torch.sum(E)
        return energy

    def compute_elastic_energy(self):
        energy = 0
        for e, tet in enumerate(self.tets):
            energy += self.compute_elem_energy(e)
        return energy

    def compute_elem_energy(self, e):
        F = self.compute_deformation_gradient(e)
        J = torch.det(F)
        if J < 0:
            return float('inf')
        C = torch.matmul(F.T, F)
        I1 = torch.trace(C)
        logJ = math.log(J)
        E = 0.5 * self.mu * (I1 - 3) - self.mu * logJ + 0.5 * self.la * logJ * logJ;
        return self.volumes[e] * E

    def compute_stress_batched(self):
        F = self.compute_deformation_gradient_batched()
        Fi = torch.inverse(F)
        Fit = torch.transpose(Fi, -1, -2)
        J = torch.det(F)
        print(J.shape, Fit.shape)
        P = self.mu * (F - Fit) + self.la * torch.log(J)[:, None, None] * Fit
        return P


    def compute_stress(self, e):
        F = self.compute_deformation_gradient(e)
        # Neo-Hookean
        Fi = torch.inverse(F)
        Fit = Fi.T
        J = torch.det(F)
        P = self.mu * (F - Fit) + self.la * math.log(J) * Fit;
        return P

    def compute_deformation_gradient_batched(self):
        xs = self.nodes[self.tets]
        d1_ = xs[:, 1] - xs[:, 0]
        d2_ = xs[:, 2] - xs[:, 0]
        d3_ = xs[:, 3] - xs[:, 0]
        D_ = torch.cat([d1_, d2_, d3_], dim=1)
        D_ = torch.transpose(D_.reshape(D_.shape[0], 3, 3), -1, -2)
        F_ = torch.bmm(D_, torch.transpose(self.torch_jacobians, -1, -2))
        return F_

    def compute_deformation_gradient(self, e):
        tet = self.tets[e]
        # get the 4 points
        x0 = self.nodes[tet[0]]
        x1 = self.nodes[tet[1]]
        x2 = self.nodes[tet[2]]
        x3 = self.nodes[tet[3]]
        # build the columns of the matrix
        d1 = x1 - x0
        d2 = x2 - x0
        d3 = x3 - x0
        # compute the deformation gradient
        D = torch.cat([d1.view(-1, *d1.shape), d2.view(-1, *d2.shape), d3.view(-1, *d3.shape)], dim=0).T
        F = torch.matmul(D, self.jacobians[e].T)
        return F

