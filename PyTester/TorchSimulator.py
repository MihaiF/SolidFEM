import torch
import numpy as np
import math

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
            Dinv = torch.inverse(D)     
            self.jacobians.append(Dinv.T)
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
        # TODO: optimize
        for i, v in enumerate(self.nodes):
            if self.inv_masses[i] != 0:
                f[i] = f[i] + self.masses[i] * self.gravity
        return f
        
    def compute_energy(self):
        energy = self.compute_elastic_energy()
        # TODO: dot product between nodes and gravity weighed by masses
        for i, v in enumerate(self.nodes):
            if self.inv_masses[i] != 0:
                energy -= self.masses[i, 0] * self.gravity[1] * self.nodes[i, 1]
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

    def compute_stress(self, e):
        F = self.compute_deformation_gradient(e)
        # Neo-Hookean
        Fi = torch.inverse(F)
        Fit = Fi.T
        J = torch.det(F)
        P = self.mu * (F - Fit) + self.la * math.log(J) * Fit;
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
        F = torch.matmul(D, self.jacobians[e].T)
        return F
    
