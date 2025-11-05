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
        reference_matrix = np.array([[ -0.331238650671255E+00,  -0.298887236553943E+00,   0.767438408492079E+00,  -0.981000000000000E+01,  -0.221359716905403E-16,   0.568545698421050E-16,  -0.368404767176652E-16,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00],
                                     [ -0.262561431782793E+01,  -0.859479839872775E+01,   0.486539857721755E+01,  -0.000000000000000E+00,  -0.105437510963931E-16,   0.105931441986074E-14,  -0.273143904057768E-15,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00],
                                     [ -0.553813632136812E+00,  -0.531457135059466E+01,  -0.110604930073433E+02,   0.000000000000000E+00,  -0.251663162816342E-15,   0.117887431061072E-14,  -0.168840984536176E-15,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00],
                                     [ -0.000000000000000E+00,  -0.000000000000000E+00,   0.100000000000000E+01,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00],
                                     [  0.188790354201620E-16,   0.304881511671334E-16,  -0.256186158055978E-16,  -0.000000000000000E+00,  -0.194900258961829E+00,  -0.146377178876203E+01,  -0.951427851116436E+01,   0.981000000000000E+01,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,  -0.000000000000000E+00],
                                     [  0.123154375350080E-14,  -0.103904882878353E-13,  -0.666859409966363E-15,  -0.000000000000000E+00,   0.265010246824567E+00,  -0.547978547561464E+02,   0.222343330761832E+02,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00],
                                     [ -0.136799131304572E-15,   0.243592681111350E-15,  -0.929689485729121E-16,  -0.000000000000000E+00,   0.625516795982231E+00,   0.113143829605573E+00,  -0.226584466905136E+01,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00],
                                     [ -0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.100000000000000E+01,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00],
                                     [  0.100000000000000E+01,   0.000000000000000E+00,  -0.000000000000000E+00,  -0.103645417057814E+01,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00],
                                     [  0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.100000000000000E+01,  -0.000000000000000E+00,  -0.000000000000000E+00,   0.103645417057814E+01,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.994614310938120E+01],
                                     [ -0.000000000000000E+00,   0.100000000000000E+01,  -0.000000000000000E+00,  -0.994614310938120E+01,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00],
                                     [ -0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.100000000000000E+01,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00,  -0.000000000000000E+00,   0.000000000000000E+00]])
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
                mask_small = np.abs(vecs_avl_sorted[idx_eig, :]) > 1e-13

                compute_vals = (np.dot(A, vecs_avl_sorted[idx_eig, :])) / vecs_avl_sorted[idx_eig]

                # eigen vectors are realy not that accurate, only use ~5e-7
                np.testing.assert_allclose(compute_vals[mask_small], vals_avl_sorted[idx_eig,], atol=5e-7, rtol=5e-7)


if __name__ == "__main__":
    unittest.main()
