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
geom_dir = os.path.join(base_dir, '..', 'geom_files')

geom_file = os.path.join(geom_dir, "aircraft.avl")
mass_file = os.path.join(geom_dir, "aircraft.mass")

class TestEigenAnalysis(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_trim_setup(self):
        vel = 10.0

        self.ovl_solver.set_trim_condition("velocity", vel)

        # internally avl will calclate what the value of CL should be to trim at the set velocity
        dens = self.ovl_solver.get_parameter("density")
        g = self.ovl_solver.get_parameter("grav.acc.")
        mass = self.ovl_solver.get_parameter("mass")
        Sref = self.ovl_solver.get_reference_data()["Sref"]
        weight = mass * g
        cl = weight / (0.5 * dens * vel**2 * Sref)

        cl_internal = self.ovl_solver.get_avl_fort_arr("CASE_R", "PARVAL")[0][5]
        np.testing.assert_allclose(cl_internal, cl, rtol=1e-14)

    def test_sys_matrix(self):
        self.ovl_solver.set_trim_condition("velocity", 10)
        self.ovl_solver.set_constraint("Elevator", "Cm", 0.00)
        self.ovl_solver.set_parameter("X cg", 0.0654)
        self.ovl_solver.execute_run()
        self.ovl_solver.execute_eigen_mode_calc()
        Amat = self.ovl_solver.get_system_matrix(in_body_axis=True)

        # fmt:off
        # this is direct output from AVL for the same case
        reference_matrix = np.array([[-0.460867888472343E+00,   0.173055728201383E-01,   0.153140775675941E+01,  -0.981000000000000E+01,   0.376671530951219E-02,  -0.491511165410724E-03,  -0.670511471689654E-02,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00],
                                     [-0.261811973218449E+01,  -0.864283107117767E+01,   0.476449167869137E+01,  -0.000000000000000E+00,  -0.348842596009870E-03,   0.455197743457654E-04,   0.620973296940547E-03,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00],
                                     [-0.522436199390934E+00,  -0.557638031175674E+01,  -0.116062556873365E+02,   0.000000000000000E+00,  -0.350642661259852E-02,   0.452655653177726E-03,   0.622747356584008E-02,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00],
                                     [-0.000000000000000E+00,  -0.000000000000000E+00,   0.100000000000000E+01,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00],
                                     [ 0.188503004913907E-16,   0.186106272094042E-16,  -0.197017879915129E-16,  -0.000000000000000E+00,  -0.216127299883635E+00,  -0.146715537953311E+01,  -0.947976161928692E+01,   0.981000000000000E+01,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,  -0.000000000000000E+00],
                                     [ 0.395701675975566E-14,   0.522521934909772E-15,   0.237798403242147E-14,  -0.000000000000000E+00,   0.267358149012866E+00,  -0.548090472873277E+02,   0.222093889208183E+02,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00],
                                     [-0.754041521309062E-17,   0.659354598794689E-15,   0.612408507949317E-15,  -0.000000000000000E+00,   0.744021530053689E+00,   0.101883215932938E+00,  -0.247683831678735E+01,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00],
                                     [-0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.100000000000000E+01,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00],
                                     [ 0.100000000000000E+01,   0.000000000000000E+00,  -0.000000000000000E+00,  -0.103933303483671E+01,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00],
                                     [ 0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.100000000000000E+01,  -0.000000000000000E+00,  -0.000000000000000E+00,   0.103933303483671E+01,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.994584269143128E+01],
                                     [-0.000000000000000E+00,   0.100000000000000E+01,  -0.000000000000000E+00,  -0.994584269143128E+01,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00],
                                     [-0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.100000000000000E+01,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00]])
        # fmt:on
        np.testing.assert_allclose(Amat, reference_matrix, atol=1e-12)

    def test_vel_sweep(self):
        # for vel in np.linspace(10, 100, 10):
        for vel in [10]:
            self.ovl_solver.set_trim_condition("velocity", vel)

            self.ovl_solver.execute_eigen_mode_calc()
            vecs_avl = self.ovl_solver.get_eigenvectors()
            vals_avl = self.ovl_solver.get_eigenvalues()
            num_eigs = len(vals_avl)

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
                mask_small = np.abs(vecs_avl_sorted[idx_eig, :]) > 5e-13

                compute_vals = (np.dot(A, vecs_avl_sorted[idx_eig, :])) / vecs_avl_sorted[idx_eig]

                # eigen vectors are realy not that accurate, only use ~5e-7
                np.testing.assert_allclose(compute_vals[mask_small], vals_avl_sorted[idx_eig,], atol=5e-7, rtol=5e-7)


if __name__ == "__main__":
    unittest.main()
