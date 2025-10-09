# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")


class TestEigenAnalysis(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file, timing=False)
        
    def test_trim_setup(self):
        
        vel = 10.0

        self.ovl_solver.set_trim_condition("velocity", vel)

        # internally avl will calclate what the value of CL should be to trim at the set velocity
        dens = self.ovl_solver.get_parameter("density")
        g = self.ovl_solver.get_parameter("grav.acc.")
        mass = self.ovl_solver.get_parameter("mass")
        Sref = self.ovl_solver.get_reference_data()["Sref"]
        weight = mass * g
        cl = weight / (0.5 * dens * vel**2*Sref)
        
        cl_internal = self.ovl_solver.get_avl_fort_arr("CASE_R", "PARVAL")[0][5]
        np.testing.assert_allclose(cl_internal, cl, rtol=1e-14)
        
    
    def test_sys_matrix(self):
        self.ovl_solver.set_trim_condition("velocity", 10)
        self.ovl_solver.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")
        self.ovl_solver.set_parameter("X cg", 0.0653909633978)
        self.ovl_solver.execute_run()
        self.ovl_solver.execute_eigen_mode_calc()
        Amat = self.ovl_solver.get_system_matrix()
        # np.set_printoptions(precision=40, suppress=False)
        # print(Amat)
        
        reference_matrix = np.array([[-1.8902384357664129e-01,  2.3002689477803040e-01, -5.7530018186675258e-01, 9.8100000000000005e+00, -1.2944346868143671e-18,  3.0836085694547899e-16, 5.7497716015888938e-17, -0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00, -0.0000000000000000e+00],
                                     [-2.4378954434082387e+00, -7.1715499059211556e+00, -6.7560778042247893e+00, 0.0000000000000000e+00, -0.0000000000000000e+00,  2.4364377334326530e-17, 7.9356312600659169e-16,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00, -0.0000000000000000e+00],
                                     [ 3.5503919975149523e-01,  4.3901675673172225e+00, -1.0884122291305133e+01, 0.0000000000000000e+00, -5.5877831966592671e-17, -6.1024867081123897e-16, 1.5780444257853037e-16,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                                     [ 0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                                     [-2.9751685342270786e-18,  1.3625678559041760e-17, -1.8268537303960584e-17,-0.0000000000000000e+00, -1.7382462626517320e-01,  7.7056552839421799e-01, 9.6479730721370238e+00,  9.8100000000000005e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00, -0.0000000000000000e+00],
                                     [-2.0183516732755471e-15, -8.5816410546466867e-15, -6.9824056173870493e-16, 0.0000000000000000e+00, -3.1605707603783356e-01, -1.4791414776917875e+02, 6.2079779415589755e+01,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                                     [ 2.5618728997302482e-16, -5.2078424884936391e-16, -8.8992027211158922e-17, 0.0000000000000000e+00, -7.5026777342273165e-01, -3.0747232776233089e+00,-2.2902836988285507e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                                     [ 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00, -1.0000000000000000e+00,-0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                                     [ 1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 8.0608278918985932e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00, -0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                                     [-0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,-0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  8.0608278918985932e-01,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  9.9674585796466637e+00],
                                     [-0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00, 9.9674585796466637e+00, -0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00, -0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00, -0.0000000000000000e+00],
                                     [ 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,-1.0000000000000000e+00, -1.0736384262997174e-27,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00]])
                                            
        
        np.testing.assert_allclose(Amat, reference_matrix, atol=1e-12)
        

    def test_vel_sweep(self):
        
        # for vel in np.linspace(10, 100, 10):
        for vel in [10]:
            self.ovl_solver.set_trim_condition("velocity", vel)
            
            
            self.ovl_solver.execute_eigen_mode_calc()
            vecs_avl = self.ovl_solver.get_eigenvectors()
            vals_avl = self.ovl_solver.get_eigenvalues()
            num_eigs =  len(vals_avl)
            
            # use the numpy eig function to get the eigenvalues and eigenvectors
            A = self.ovl_solver.get_system_matrix()
            vals_np, vecs_np = np.linalg.eig(A)
            vecs_np = vecs_np.T
            
            # because the eigenvalues are not in a consistent order, we need to sort them
            # we will sort them by the largest magnitude
            # we will also sort the eigenvectors to match the eigenvalues
            idx_np = np.argsort(np.abs(vals_np))[::-1]
            vals_np_sorted = vals_np[idx_np]
            vecs_np_sorted = vecs_np[idx_np, :]

            idx_avl = np.argsort(np.abs(vals_avl))[::-1]
            vals_avl_sorted = vals_avl[idx_avl]
            vecs_avl_sorted = vecs_avl[idx_avl, :]
            
            
            
            np.testing.assert_allclose(vals_avl_sorted, vals_np_sorted[:num_eigs], rtol=5e-14)


            # the eigenvecs appear to be poorly conditioned. 
            # OptVL, AVL, and numpy are all slightly different.
            # The almost zero values of the A matrix input to the eigen solver can be different to 3E-015 vs 2.E-015
            # I traced back the differences between OptVL and AVL through the Eigsol and appears the pivioting and scalling at one point is differenct due to rounding errors
            # Just verify that each vector is an acutal eigenvector
            for idx_eig in range(num_eigs):
                
                # only check the values of the eigen vectors greater thann 1e-14
                mask_small = np.abs(vecs_avl_sorted[idx_eig, :]) > 1e-13
                
                compute_vals = (np.dot(A,vecs_avl_sorted[idx_eig, :]))/vecs_avl_sorted[idx_eig]
                
                # eigen vectors are realy not that accurate, only use ~5e-7
                np.testing.assert_allclose(compute_vals[mask_small], vals_avl_sorted[idx_eig,], atol=5e-7,rtol=5e-7)
                
if __name__ == "__main__":
    unittest.main()
