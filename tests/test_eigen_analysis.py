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


class TestEigenAnalysisSweep(unittest.TestCase):
    def setUp(self):
        self.avl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file, timing=False)

    def test_vel_sweep(self):
        
        # for vel in np.linspace(10, 100, 10):
        for vel in [10]:
            self.avl_solver.set_parameter("velocity", vel)
            dens = self.avl_solver.get_parameter("density")
            g = self.avl_solver.get_parameter("grav.acc.")
            mass = self.avl_solver.get_parameter("mass")
            weight = mass * g
            cl = weight / (0.5 * dens * vel**2)
            self.avl_solver.set_trim_condition("CL", cl)
            
            
            self.avl_solver.execute_eigen_mode_calc()
            vecs_avl = self.avl_solver.get_eigenvectors()
            vals_avl = self.avl_solver.get_eigenvalues()
            num_eigs =  len(vals_avl)
            
            # use the numpy eig function to get the eigenvalues and eigenvectors
            A = self.avl_solver.get_system_matrix()
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
