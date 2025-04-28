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
import numpy as np
import sys

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestImport(unittest.TestCase):
    # TODO: add test for expected input output errors
    def test_instances(self):

        from optvl import OVLSolver

        avl_solver1 = OVLSolver(geo_file=geom_file)

        avl_solver2 = OVLSolver(geo_file=geom_file)

        assert avl_solver1.avl is not avl_solver2.avl

        avl_solver1.set_avl_fort_arr("CASE_R", "ALFA", 1.1)
        avl_solver2.set_avl_fort_arr("CASE_R", "ALFA", 2.0)

        assert avl_solver1.get_avl_fort_arr("CASE_R", "ALFA") != avl_solver2.get_avl_fort_arr("CASE_R", "ALFA")


if __name__ == "__main__":
    unittest.main()
