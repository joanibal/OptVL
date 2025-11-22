# =============================================================================
# Extension modules
# =============================================================================

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, '..', 'geom_files')

geom_file = os.path.join(geom_dir, "aircraft.avl")
mass_file = os.path.join(geom_dir, "aircraft.mass")


class TestImport(unittest.TestCase):
    # TODO: add test for expected input output errors
    def test_instances(self):
        from optvl import OVLSolver

        ovl_solver1 = OVLSolver(geo_file=geom_file)

        ovl_solver2 = OVLSolver(geo_file=geom_file)

        assert ovl_solver1.avl is not ovl_solver2.avl

        ovl_solver1.set_avl_fort_arr("CASE_R", "ALFA", 1.1)
        ovl_solver2.set_avl_fort_arr("CASE_R", "ALFA", 2.0)

        assert ovl_solver1.get_avl_fort_arr("CASE_R", "ALFA") != ovl_solver2.get_avl_fort_arr("CASE_R", "ALFA")


if __name__ == "__main__":
    unittest.main()
