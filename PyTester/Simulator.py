import numpy as np
import math
import scipy
from scipy.optimize import minimize

class Simulator:

    gravity = np.array([0, -9.8, 0])
    density = 1070

    def __init__(self,
                nodes,
                tets,
                fixed_nodes,
                config
        ):
        """
        @param nodes: Array of node positions
        @param tets: Array of tetrahedra indices
        @param fixed_nodes: Array of indices indicating fixed nodes
        @param config: Dictionary containing configuration
        """        
        self.substeps = 1
        if "substeps" in config:
            self.substeps = config["substeps"]
        self.dt = 0.016 / self.substeps
            
        # init members
        self.nodes = nodes.astype('double')
        numTets = int(tets.shape[0] / 4)
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
        self.young = config["young"]
        self.poisson = config["poisson"]
        self.mu = 0.5 * self.young / (1 + self.poisson)
        self.la = self.young * self.poisson / (1 + self.poisson) / (1 - 2 * self.poisson)        
        
    def step(self):
        for iter in range(0, self.substeps):
            self.step_explicit()
        return self.nodes
        
    def step_explicit(self):  
        # compute forces
        f = self.compute_forces()
        # Semi-implicit (symplectic) Euler integration        
        self.vel[self.free] = self.vel[self.free] + self.dt * (self.gravity + self.inv_masses[self.free] * f[self.free])
        self.nodes[self.free] = self.nodes[self.free] + self.dt * self.vel[self.free]
        return self.nodes
    
    def step_static(self):
        return self.solve_scipy()

    def solve_scipy(self, solver='Newton-CG', verbose=False):
        x0 = self.nodes[self.free].flatten()
        sol = minimize(self.energy_objective, x0, method=solver, jac=self.jacobian)
        if verbose:
            print(sol.message)
        self.nodes[self.free] = np.reshape(sol.x, (int(sol.x.shape[0] / 3), 3))
        #print('E: ', self.compute_energy())
        return self.nodes
        
    def energy_objective(self, x):
        self.nodes[self.free] = np.reshape(x, (int(x.shape[0] / 3), 3))
        energy = self.compute_energy()
        #print('energy ', energy)
        return energy
    
    def jacobian(self, x):
        self.nodes[self.free] = np.reshape(x, (int(x.shape[0] / 3), 3))
        g = -self.compute_gradients()
        return g[self.free].flatten()
        
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
                print('Not converging after {0} iterations'.format(iter+1))
                break
        print(energy)
        return self.nodes

    def solve_steep_desc(self, num_iters=100):
        for iter in range(0, num_iters):
            print('E: ', self.compute_energy())
            r = self.compute_gradients()
            dx = np.zeros(r.shape)
            dx[self.free] = r[self.free]
            df = -self.compute_force_diff(dx)
            delta = np.dot(r[self.free].flatten(), r[self.free].flatten())
            alpha = delta / np.dot(r[self.free].flatten(), df[self.free].flatten())
            self.nodes[self.free] = self.nodes[self.free] + alpha * r[self.free]
        return self.nodes

    def solve_conj_grad(self, num_iters=10, ls_iter=2, ls_tol=1e-10):
        r = self.compute_gradients()
        d = r
        delta = np.dot(r[self.free].flatten(), r[self.free].flatten())
        delta0 = delta
        for iter in range(0, num_iters):
            deltaD = np.dot(d[self.free].flatten(), d[self.free].flatten())
            for j in range(0, ls_iter):
                dx = np.zeros(d.shape)
                dx[self.free] = d[self.free]
                df = -self.compute_force_diff(dx)
                a = np.dot(r[self.free].flatten(), d[self.free].flatten())
                b = np.dot(d[self.free].flatten(), df[self.free].flatten())
                alpha = a / b
                self.nodes[self.free] = self.nodes[self.free] + alpha * d[self.free]
                r = self.compute_gradients()
                #print('tol: ', math.sqrt(alpha * alpha * deltaD))
                if alpha * alpha * deltaD < ls_tol * ls_tol:
                    #print('line search stopped at iteration ', j + 1)
                    break            
            delta1 = np.dot(r[self.free].flatten(), r[self.free].flatten())
            beta = delta1 / delta
            delta = delta1
            d = r + beta * d;
            #print('E: ', self.compute_energy(), delta)
            if delta < 1e-5:
                print('Converged at ', iter + 1)
                return
        print('Not converged')

    def estimate_params(self, target):
        self.target = target
        x0 = [self.mu, self.la]
        sol = minimize(self.inverse_objective, x0, method='Powell', bounds=((0, None), (0, None)))
        print(sol.message)
        self.mu = sol.x[0]
        self.la = sol.x[1]
        print(self.mu, self.la)        

    def inverse_objective(self, x):
        self.mu = x[0]
        self.la = x[1]
        print(self.mu, self.la)
        print('Solving...')
        #self.solve_scipy(solver='BFGS')
        self.solve_conj_grad(400)
        delta = (self.nodes - self.target).flatten()
        error = np.dot(delta, delta)
        print('err {:E}'.format(error))
        return error

    # returns the elastic forces
    def compute_forces(self):
        f = np.zeros(self.nodes.shape)
        for e, tet in enumerate(self.tets):
            P = self.compute_stress(e)
            H = -self.volumes[e] * np.matmul(P, self.jacobians[e]) 
            for j in range(1, 4):
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
            if J < 0:
                return float('inf')
            C = np.matmul(np.transpose(F), F)
            I1 = np.trace(C)
            #I3 = np.linalg.det(C)
            #J = math.sqrt(I3)
            logJ = math.log(J)
            E = 0.5 * self.mu * (I1 - 3) - self.mu * logJ + 0.5 * self.la * logJ * logJ
            energy += self.volumes[e] * E
        return energy
            
    def compute_stress(self, e):
        F = self.compute_deformation_gradient(e)
        # Neo-Hookean
        Fi = np.linalg.inv(F)
        Fit = np.transpose(Fi)
        J = np.linalg.det(F)
        P = self.mu * (F - Fit) + self.la * math.log(J) * Fit
        return P
        
    def compute_force_diff(self, dx):
        df = np.zeros(self.nodes.shape)
        dx_local = np.empty((4, 3))
        for e, tet in enumerate(self.tets):
            dx_local[0] = dx[tet[0]]
            dx_local[1] = dx[tet[1]]
            dx_local[2] = dx[tet[2]]
            dx_local[3] = dx[tet[3]]
            df_local = self.compute_local_force_diff(e, dx_local)
            df[tet[0]] += df_local[0]
            df[tet[1]] += df_local[1]
            df[tet[2]] += df_local[2]
            df[tet[3]] += df_local[3]
        return df

    def compute_local_force_diff(self, e, dx):
        tet = self.tets[e]
        # get the 4 differential points
        x0 = dx[0]
        x1 = dx[1]
        x2 = dx[2]
        x3 = dx[3]
        # build the columns of the matrix
        d1 = x1 - x0
        d2 = x2 - x0
        d3 = x3 - x0
        # compute the deformation gradient differential
        D = np.column_stack((d1, d2, d3))
        dF = np.matmul(D, np.transpose(self.jacobians[e]))
        # the deformation gradient and related quantities
        F = self.compute_deformation_gradient(e)
        Fi = np.linalg.inv(F)
        Fit = np.transpose(Fi)
        J = np.linalg.det(F)
        # the stress differential
        A = np.matmul(Fi, dF)
        B = np.matmul(Fit, np.transpose(A))        
        dP = self.mu * dF + (self.mu - self.la * math.log(J)) * B + self.la * np.trace(A) * Fit
        dH = -self.volumes[e] * np.matmul(dP, self.jacobians[e])
        df = np.empty((4, 3))
        df[1] = dH[:,0]
        df[2] = dH[:,1]
        df[3] = dH[:,2]
        df[0] = -(df[1] + df[2] + df[3])
        return df
        
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