import numpy as np
import math

class Simulator:

    dt = 0.016
    gravity = np.array([0, -9.8, 0])
    young = 66000
    poisson = 0.4
    density = 1070
    volumes = []
    jacobians = []

    def __init__(self,
                nodes,
                tets,
                particle_mass
        ):
        """
        @param nodes: Array of node positions
        @param tets: Array of tetrahedra indices
        """        
        # init members
        self.nodes = nodes
        numTets = math.ceil(tets.shape[0] / 4)
        self.tets = tets.reshape((numTets, 4))
        self.masses = np.zeros((nodes.shape[0], 1), dtype=np.double)
        self.vel = np.zeros(nodes.shape, dtype=np.double)
        # compute volumes and masses
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
        for i, m in enumerate(self.masses):            
            if m != 0 and particle_mass[i] != 0:
                self.inv_masses[i] = 1.0 / m
        # init params
        self.mu = 0.5 * self.young / (1 + self.poisson)
        self.la = self.young * self.poisson / (1 + self.poisson) / (1 - 2 * self.poisson)        
        
    def step(self):
        return self.step_explicit()
        
    def step_explicit(self):  
        # compute forces
        f = self.compute_gradients()
        # Semi-implicit (symplectic) Euler integration        
        for i, v in enumerate(self.nodes):
            if self.inv_masses[i] != 0:
                self.vel[i] = self.vel[i] + self.dt * (self.gravity + self.inv_masses[i] * f[i])
                self.nodes[i] = self.nodes[i] + self.dt * self.vel[i]
        return self.nodes
        
    def compute_gradients(self):
        f = np.zeros(self.nodes.shape)
        for e, tet in enumerate(self.tets):
            P = self.compute_stress(e)
            H = -self.volumes[e] * np.matmul(P, self.jacobians[e]) 
            for j in range(1,4):
                force = H[:, j - 1]
                f[tet[j]] += force
                f[tet[0]] -= force
        return f
        
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