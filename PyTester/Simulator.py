import numpy as np
import math

class Simulator:

    dt = 0.016
    gravity = np.array([0, -9.8, 0])
    young = 66000
    poisson = 0.4
    density = 1070

    def __init__(self,
                nodes,
                tets,
                fixed_nodes
        ):
        """
        @param nodes: Array of node positions
        @param tets: Array of tetrahedra indices
        @param fixed_nodes: Array of indices indicating fixed nodes
        """        
        # init members
        self.nodes = nodes
        numTets = math.ceil(tets.shape[0] / 4)
        self.tets = tets.reshape((numTets, 4))
        self.masses = np.zeros((nodes.shape[0], 1), dtype=np.double)
        self.vel = np.zeros(nodes.shape, dtype=np.double)
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
            D = np.column_stack((d1, d2, d3))
            vol = np.linalg.det(D) / 6.0
            self.volumes.append(vol)
            Dinv = np.linalg.inv(D)            
            self.jacobians.append(np.transpose(Dinv))
            # compute the element mass
            mass = 0.25 * vol * self.density
            self.masses[tet[0]] += mass;
            self.masses[tet[1]] += mass;
            self.masses[tet[2]] += mass;
            self.masses[tet[3]] += mass;
        # compute inverse masses
        self.inv_masses = np.zeros(self.masses.shape, dtype=np.double)
        self.free = np.setdiff1d(range(0, nodes.shape[0]), fixed_nodes)
        self.inv_masses[self.free] = 1.0 / self.masses[self.free]
        # init params
        self.mu = 0.5 * self.young / (1 + self.poisson)
        self.la = self.young * self.poisson / (1 + self.poisson) / (1 - 2 * self.poisson)        
        
    def step(self):
        return self.step_explicit()
        
    def step_explicit(self):  
        # compute forces
        f = self.compute_forces()
        # Semi-implicit (symplectic) Euler integration        
        self.vel[self.free] = self.vel[self.free] + self.dt * (self.gravity + self.inv_masses[self.free] * f[self.free])
        self.nodes[self.free] = self.nodes[self.free] + self.dt * self.vel[self.free]
        return self.nodes
    
    def step_static(self):
        # gradient descent
        numIters = 100
        alpha = 1e-3
        energy = self.compute_energy()
        for iter in range(0, numIters):
            f = self.compute_gradients()            
            self.nodes[self.free] = self.nodes[self.free] + alpha * f[self.free]
            old_energy = energy;
            energy = self.compute_energy()
            rel_err = abs((energy - old_energy) / energy)
            if rel_err < 1e-15:
                print('Converged after {0} iterations'.format(iter+1))
                break
        return self.nodes

    # returns the elstic forces
    def compute_forces(self):
        f = np.zeros(self.nodes.shape)
        for e, tet in enumerate(self.tets):
            P = self.compute_stress(e)
            H = -self.volumes[e] * np.matmul(P, self.jacobians[e]) 
            for j in range(1,4):
                force = H[:, j - 1]
                f[tet[j]] += force
                f[tet[0]] -= force
        return f
        
    # returns the negative gradients
    def compute_gradients(self):
        f = self.compute_forces()
        for i, v in enumerate(self.nodes):
            if self.inv_masses[i] != 0:
                f[i] = f[i] + self.masses[i] * self.gravity
        return f
        
    def compute_energy(self):
        energy = self.compute_elastic_energy()
        # TODO: dot product between nodes and gravity weighed by masses
        for i, v in enumerate(self.nodes):
            if self.inv_masses[i] != 0:
                energy -= self.masses[i] * self.gravity[1] * self.nodes[i, 1]
        return energy
        
    def compute_elastic_energy(self):
        energy = 0
        for e, tet in enumerate(self.tets):
            F = self.compute_deformation_gradient(e)
            J = np.linalg.det(F)
            C = np.matmul(np.transpose(F), F)
            I1 = np.trace(C)
            logJ = math.log(J)
            E = 0.5 * self.mu * (I1 - 3) - self.mu * logJ + 0.5 * self.la * logJ * logJ;
            energy += self.volumes[e] * E
        return energy
            
    def compute_stress(self, e):
        F = self.compute_deformation_gradient(e)
        # Neo-Hookean
        Fi = np.linalg.inv(F)
        Fit = np.transpose(Fi)
        J = np.linalg.det(F)
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
        D = np.column_stack((d1, d2, d3))        
        F = np.matmul(D, np.transpose(self.jacobians[e]))
        return F