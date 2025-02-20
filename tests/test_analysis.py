# =============================================================================
# Extension modules
# =============================================================================
from optvl import AVLSolver

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


class TestAnalysisSweep(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file, timing=False)

    def test_constrained_alpha_sweep(self):
        self.avl_solver.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")
        self.avl_solver.set_constraint("Rudder", 0.00, con_var="Cn yaw moment")

        alpha_array = np.arange(0, 10)
        cl_ref_arr = np.array(
            [
                1.2349061238204873,
                1.3282779560514677,
                1.42113603863798,
                1.5134418218893897,
                1.6051572516385355,
                1.6962448080698307,
                1.7866675434128751,
                1.8763891184181036,
                1.9653738375342211,
                2.053586682710641,
            ]
        )
        cd_ref_arr = np.array(
            [
                0.030753634421301752,
                0.033676747268002544,
                0.03679211508827116,
                0.040092804096450746,
                0.04357129507051221,
                0.047219503496354835,
                0.05102880128527856,
                0.0549900400052339,
                0.05909357556230493,
                0.06332929426488881,
            ]
        )
        for idx_alpha, alpha in enumerate(alpha_array):
            self.avl_solver.set_constraint("alpha", alpha)

            self.avl_solver.execute_run()
            run_data = self.avl_solver.get_total_forces()

            np.testing.assert_allclose(
                cl_ref_arr[idx_alpha],
                run_data["CL"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                cd_ref_arr[idx_alpha],
                run_data["CD"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                run_data["CM"],
                0.0,
                atol=1e-8,
            )

    def test_constrained_cl_sweep(self):
        self.avl_solver.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")
        self.avl_solver.set_constraint("Rudder", 0.00, con_var="Cn yaw moment")

        cd_ref_arr = np.array(
            [
                0.01654814255244833,
                0.018124778345018383,
                0.01994896285331091,
                0.02202067604738557,
                0.024339528639367926,
                0.026904749658315106,
                0.029715171166770832,
                0.03276920986323026,
                0.03606484526169707,
                0.03959959407831794,
                0.043370480383046583,
            ]
        )
        cl_arr = np.arange(0.6, 1.7, 0.1)
        for idx_cl, cl in enumerate(cl_arr):
            self.avl_solver.set_trim_condition("CL", cl)
            # the tight tolerance here helps catch small issues with the newton solver
            self.avl_solver.execute_run(tol=1e-12)
            run_data = self.avl_solver.get_total_forces()

            np.testing.assert_allclose(
                cl,
                run_data["CL"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                cd_ref_arr[idx_cl],
                run_data["CD"],
                rtol=1e-8,
            )
            np.testing.assert_allclose(
                run_data["CM"],
                0.0,
                atol=1e-8,
            )


class TestBodyAnalysis(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file="supra.avl",debug=False)
    
    def test_coefs(self):
        self.avl_solver.set_constraint("alpha", 5.00)
        self.avl_solver.execute_run()
        coef_data = self.avl_solver.get_total_forces()

        # the values are wonky here because of an unrealistic CDCL curve
        np.testing.assert_allclose(coef_data["CL"], 0.636031170179549, rtol=1e-8)
        np.testing.assert_allclose(coef_data["CD"], 3.6953247032454204, rtol=1e-8)
        np.testing.assert_allclose(coef_data["CM"], -0.5736410313952236, rtol=1e-8)

class TestHingeMom(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_con_surf_mom(self):
        self.avl_solver.set_constraint("Elevator", 10.00)
        self.avl_solver.set_constraint("Rudder", 0.00)
        self.avl_solver.execute_run()

        mom_data = self.avl_solver.get_hinge_moments()

        np.testing.assert_allclose(mom_data["Elevator"], -0.04381216304, rtol=1e-8)
        np.testing.assert_allclose(mom_data["Rudder"], 0.0, atol=1e-8)

        self.avl_solver.set_constraint("Elevator", 0.00)
        self.avl_solver.set_constraint("Rudder", 10.00)
        self.avl_solver.execute_run()
        mom_data = self.avl_solver.get_hinge_moments()

        np.testing.assert_allclose(mom_data["Rudder"], -1.0957068091302286e-3, rtol=1e-8)
        np.testing.assert_allclose(mom_data["Elevator"], -1.065702741142345e-2, rtol=1e-8)


class TestCaseDerivs(unittest.TestCase):
    def setUp(self) -> None:
        # self.avl_solver = AVLSolver(geo_file=geom_file)
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")

    def test_coefs_wrt_con_surfs(self):
        self.avl_solver.set_constraint("alpha", 45.00)
        self.avl_solver.execute_run()
        run_data = self.avl_solver.get_total_forces()
        coef_derivs = self.avl_solver.get_control_stab_derivs()
        # TODO: test againast values from AVL


class TestVariableSetting(unittest.TestCase):
    def setUp(self) -> None:
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_alpha_set(self):
        """
        Test that setting the alpha works
        """

        self.avl_solver.set_constraint("alpha", 10.00)
        self.avl_solver.execute_run()
        alfa_new = self.avl_solver.get_parameter("alpha")
        np.testing.assert_allclose(alfa_new, 10.00, rtol=1e-15)

    def test_con_surf_set(self):
        """
        Test that setting the control surface works
        """

        self.avl_solver.set_constraint("Elevator", 10.00)
        self.avl_solver.set_constraint("alpha", 10.00)
        self.avl_solver.execute_run()
        alfa_new = self.avl_solver.get_parameter("alpha")
        np.testing.assert_allclose(alfa_new, 10.00, rtol=1e-15)
        def_dict = self.avl_solver.get_control_deflections()
        np.testing.assert_allclose(def_dict["Elevator"], 10.00, rtol=1e-15)


if __name__ == "__main__":
    unittest.main()
