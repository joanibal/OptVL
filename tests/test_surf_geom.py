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
geom_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestGeom(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file=geom_file)

    def test_surface_params(self):
        reference_data = {
            "Wing": {
                "nchordwise": 7,
                "cspace": 1.0,
                "nspan": 20,
                "sspace": -2.0,
                "yduplicate": 0.0,
                "scale": np.array([1.1, 1.2, 1.3]),
                "translate": np.array([0.1, 0.2, 0.3]),
                "angle": 1.23,
                "nspans": np.array([5, 4, 3, 2, 1]),
                "sspaces": np.array([-1.0, 0.0, 1.0, 2.0, 3.0]),
                "aincs": np.array([0.5, 0.4, 0.3, 0.2, 0.1]),
                "chords": np.array([0.5, 0.4, 0.3, 0.2, 0.1]),
                "xyzles": np.array([[0, 0, 0], [0.1, 1.0, 0.01], [0.2, 2.0, 0.02], [0.3, 3.0, 0.03], [0.4, 4.0, 0.04]]),
            },
        }

        data = self.avl_solver.get_surface_params(include_geom=True, include_panneling=True, include_con_surf=True)
        
        from pprint import pprint
        

        for surf in reference_data:
            for key in reference_data[surf]:
                np.testing.assert_allclose(
                    data[surf][key],
                    reference_data[surf][key],
                    atol=1e-8,
                    err_msg=f"Surface `{surf}` key `{key}` does not match reference data",
                )

        self.avl_solver.set_constraint("alpha", 6.00)
        self.avl_solver.set_constraint("beta", 2.00)
        self.avl_solver.execute_run()
        
        assert self.avl_solver.get_num_surfaces() == 5
        assert self.avl_solver.get_num_strips() == 90
        assert self.avl_solver.get_mesh_size() == 780

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
        
        coefs = self.avl_solver.get_total_forces()
        np.testing.assert_allclose(
            coefs["CL"],
            5.407351081559913,
            rtol=1e-8,
        )

        self.avl_solver.set_surface_params(data)
        

        assert self.avl_solver.get_num_surfaces() == 5
        assert self.avl_solver.get_num_strips() == 90
        assert self.avl_solver.get_mesh_size() == 780

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
        
        coefs = self.avl_solver.get_total_forces()
        np.testing.assert_allclose(
            coefs["CL"],
            5.407351081559913,
            rtol=1e-8,
        )

    def test_surface_mirroring(self):
        
        # Take the one wing and streach out the tip
        new_data = {
            "Wing": {
                "xyzles": np.array([[0, 0, 0], [0.1, 1.0, 0.01], [0.2, 2.0, 0.02], [0.3, 3.0, 0.03], [0.4, 10.0, 0.04]]),
            },
        }
        
        self.avl_solver.set_constraint("alpha", 10.00)
        self.avl_solver.set_surface_params(new_data)
        
        self.avl_solver.execute_run()
        
        run_data = self.avl_solver.get_total_forces()
        
        # if only one wing was updated then will have unbalanced yaw and roll moments
        np.testing.assert_allclose(
            run_data["CR SA"],
            0.0,
            atol=1e-12
        )
        
        np.testing.assert_allclose(
            run_data["CN SA"],
            0.0,
            atol=1e-12
        )
        
        updated_data = self.avl_solver.get_surface_params(include_geom=True)
        
        np.testing.assert_allclose(
            updated_data["Wing"]["xyzles"],
            new_data["Wing"]["xyzles"],
            atol=1e-16
        )


if __name__ == "__main__":
    unittest.main()
