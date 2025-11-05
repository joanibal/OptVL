# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

# =============================================================================
# Standard Python Modules
# =============================================================================
import os
import psutil

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, '..', 'geom_files')

geom_file = os.path.join(geom_dir, "aircraft.avl")
mass_file = os.path.join(geom_dir, "aircraft.mass")
geom_mod_file = os.path.join(geom_dir, "aircraft_mod.avl")
geom_output_file = os.path.join(geom_dir, "aircraft_out.avl")

supra_geom_file = os.path.join(geom_dir, 'supra.avl')
rect_geom_file = os.path.join(geom_dir, "rect.avl")
rect_geom_output_file = os.path.join(geom_dir, "rect_out.avl")

# TODO: add test for expected input output errors

class TestInput(unittest.TestCase):
    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id():80} Memory usage: {mb_memory:.2f} MB")

    def test_read_geom(self):
        ovl_solver = OVLSolver(geo_file=geom_file)
        assert ovl_solver.get_num_surfaces() == 5
        assert ovl_solver.get_num_strips() == 90
        assert ovl_solver.get_mesh_size() == 780

    def test_read_geom_and_mass(self):
        ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        assert ovl_solver.get_avl_fort_arr("CASE_L", "LMASS")


class TestOutput(unittest.TestCase):
    
    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id():80} Memory usage: {mb_memory:.2f} MB")

    
    def test_write_geom(self):
        """check that the file written by OptVL is the same as the original file"""
        ovl_solver = OVLSolver(geo_file=supra_geom_file)
        ovl_solver.write_geom_file(geom_output_file)
        baseline_data = ovl_solver.get_surface_params()
        baseline_data_body = ovl_solver.get_body_params()

        del ovl_solver
        ovl_solver = OVLSolver(geo_file=geom_output_file)
        new_data = ovl_solver.get_surface_params()
        new_data_body = ovl_solver.get_body_params()

        for surf in baseline_data:
            for key in baseline_data[surf]:
                data = new_data[surf][key]
                # check if it is a list of strings
                if isinstance(data, list) and isinstance(data[0], str):
                    for a, b in zip(data, baseline_data[surf][key]):
                        assert a == b
                else:
                    np.testing.assert_allclose(
                        new_data[surf][key],
                        baseline_data[surf][key],
                        atol=1e-8,
                        err_msg=f"Surface `{surf}` key `{key}` does not match reference data",
                    )
    
        for body in baseline_data_body:
            for key in baseline_data_body[body]:
                data = new_data_body[body][key]
                # check if it is a list of strings
                if isinstance(data, str):
                    assert new_data_body[body][key] == baseline_data_body[body][key]
                else:
                    np.testing.assert_allclose(
                        new_data_body[body][key],
                        baseline_data_body[body][key],
                        atol=1e-8,
                        err_msg=f"bodyace `{body}` key `{key}` does not match reference data",
                    )
    
    def test_write_panneling_params(self):
        # test that the surface is output correctly when only section or surface
        # panneling is given
        ovl_solver = OVLSolver(geo_file=rect_geom_file)
        ovl_solver.write_geom_file(rect_geom_output_file)   
        baseline_data = ovl_solver.get_surface_params(include_paneling=True, include_geom=False)

        assert baseline_data['Wing']['use surface spacing'] == True
        
        del ovl_solver
        ovl_solver = OVLSolver(geo_file=rect_geom_output_file)
        new_data = baseline_data = ovl_solver.get_surface_params()

        for surf in baseline_data:
            for key in baseline_data[surf]:
                data = new_data[surf][key]
                # check if it is a list of strings
                if isinstance(data, list) and isinstance(data[0], str):
                    for a, b in zip(data, baseline_data[surf][key]):
                        assert a == b
                else:
                    np.testing.assert_allclose(
                        new_data[surf][key],
                        baseline_data[surf][key],
                        atol=1e-8,
                        err_msg=f"Surface `{surf}` key `{key}` does not match reference data",
                    )


class TestFortranLevelAPI(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_get_scalar(self):
        avl_version = 3.52
        version = self.ovl_solver.get_avl_fort_arr("CASE_R", "VERSION")
        self.assertEqual(version, avl_version)

        # test that this works with lower case
        version = self.ovl_solver.get_avl_fort_arr("case_r", "version")
        self.assertEqual(version, avl_version)

    def test_get_array(self):
        chords = self.ovl_solver.get_avl_fort_arr("SURF_GEOM_R", "CHORDS")

        self.assertEqual(chords.shape, (100, 301))
        np.testing.assert_array_equal(chords[0, :5], np.array([0.45, 0.45, 0.4, 0.3, 0.2]))


if __name__ == "__main__":
    unittest.main()
