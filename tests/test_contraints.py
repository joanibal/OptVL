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


class TestConstraints(unittest.TestCase):
    # TODO: add reference values for comparison by running avl binnary by hand or with wrappper

    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)

    # TODO: add test for roll, pitch, and yaw rate
    # def test_rates(self):

    def test_angles(self):
        self.avl_solver.set_constraint("alpha", 6.00)
        self.avl_solver.set_constraint("beta", 2.00)
        self.avl_solver.execute_run()

        np.testing.assert_allclose(
            self.avl_solver.get_constraint("alpha"),
            6.0,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            self.avl_solver.get_constraint("beta"),
            2.0,
            rtol=1e-8,
        )
        run_data = self.avl_solver.get_total_forces()
        np.testing.assert_allclose(
            run_data["CL"],
            1.83005581269135,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CD"],
            0.054208182783309425,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CM"],
            -0.19528584013607497,
            rtol=1e-8,
        )

    def test_control_surfaces(self):
        self.avl_solver.set_constraint("D1", 10.00)
        self.avl_solver.set_constraint("D2", 5.00)
        self.avl_solver.execute_run()
        run_data = self.avl_solver.get_total_forces()
        
        np.testing.assert_allclose(
            run_data["CL"],
            1.0106168310619361,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CD"],
            0.03876013948341263,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CM"],
            0.959676732257107,
            rtol=1e-8,
        )

    def test_control_surfaces_names(self):
        """test that the control surface names are can be used as well"""
        self.avl_solver.set_constraint("Elevator", 10.00)
        self.avl_solver.set_constraint("Rudder", 5.00)
        self.avl_solver.execute_run()
        run_data = self.avl_solver.get_total_forces()
        np.testing.assert_allclose(
            run_data["CL"],
            1.0106168310619361,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CD"],
            0.03876013948341263,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CM"],
            0.959676732257107,
            rtol=1e-8,
        )

    def test_trim(self):
        # TODO: parametrize for the other options
        self.avl_solver.set_trim_condition("CL", 1.0)
        self.avl_solver.execute_run()
        np.testing.assert_allclose(
            self.avl_solver.get_parameter("alpha"),
            -1.8615972119506075,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            self.avl_solver.get_parameter("beta"),
            0.0,
            rtol=1e-8,
        )

        run_data = self.avl_solver.get_total_forces()

        np.testing.assert_allclose(
            run_data["CL"],
            1.0,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CD"],
            0.025388893842317614,
            rtol=1e-8,
        )
        np.testing.assert_allclose(
            run_data["CM"],
            0.2558492940049343,
            rtol=1e-8,
        )


if __name__ == "__main__":
    unittest.main()
