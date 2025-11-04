# Extension modules
from optvl import OVLSolver
from ovl_avl_comparison_utils import run_comparison

# Standard Python modules
import os
import json

# External Python modules
import unittest
import numpy as np

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, '..', 'geom_files')

geom_file = os.path.join(geom_dir, "aircraft.avl")
mod_geom_file = os.path.join(geom_dir, "aircraft_mod.avl")
mass_file = os.path.join(geom_dir, "aircraft.mass")

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
        "xles" : np.array([0, 0.1, 0.2, 0.3, 0.4]),
        "yles" : np.array([0, 1.0, 2.0, 3.0, 4.0]),
        "zles" : np.array([0, 0.01, 0.02, 0.03, 0.04])
    },
}


class TestGeom(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=mod_geom_file, mass_file=mass_file)

    def test_get_surface_params(self):

        data = self.ovl_solver.get_surface_params(include_geom=True, include_paneling=True, include_con_surf=True)
        
        for surf in reference_data:
            for key in reference_data[surf]:
                np.testing.assert_allclose(
                    data[surf][key],
                    reference_data[surf][key],
                    rtol=1e-15,
                    err_msg=f"Surface `{surf}` key `{key}` does not match reference data",
                )
    
    def test_set_surface_params(self):
        """ set the modifiied geom to have the same parameters as the unmodified and check that the analysis after is the same """
        ovl_orig = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        surf_data_orig = ovl_orig.get_surface_params()
        
        self.ovl_solver.set_surface_params(surf_data_orig)
        
        with open("avl_analysis_references/constrained_aircraft.json") as f:
            ref_data_cases = json.load(f)
            run_comparison(self.ovl_solver, ref_data_cases, rtol=1e-15, atol=1e-14)

    def test_surface_mirroring(self):
        
        # Take the one wing and streach out the tip
        new_data = {
            "Wing": {
                "yles" : np.array([0, 1.0, 2.0, 3.0, 10.0]),
            },
        }
        
        self.ovl_solver.set_variable("alpha", 10.00)
        self.ovl_solver.set_surface_params(new_data)
        
        self.ovl_solver.execute_run()
        
        run_data = self.ovl_solver.get_total_forces()
        
        # if only one wing was updated then will have unbalanced yaw and roll moments
        np.testing.assert_allclose(
            run_data["Cl"],
            0.0,
            atol=1e-12
        )
        
        np.testing.assert_allclose(
            run_data["Cn"],
            0.0,
            atol=1e-12
        )
        
        updated_data = self.ovl_solver.get_surface_params(include_geom=True)
        
        np.testing.assert_allclose(
            updated_data["Wing"]["yles"],
            new_data["Wing"]["yles"],
            atol=5e-16
        )


if __name__ == "__main__":
    unittest.main()
