# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

# =============================================================================
# Standard Python Modules
# =============================================================================
import os
import copy

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, 'geom_files')
geom_file = os.path.join(geom_dir, "aircraft.avl")
mass_file = os.path.join(geom_dir, "aircraft.mass")

class TestParameterAPI(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.params_baseline = {
            # "bank": 0.0,
            # "elevation": 0,
            "Mach": 0.1,
            "velocity": 0.0,
            "density": 1.225,
            "grav.acc.": 9.81,
            "mass": 5.63,
            "Ixx": 0.196,
            "Iyy": 1.54,
            "Izz": 1.74,
            "X cg": 0.0654,
            "Y cg": 0,
            "Z cg": 0,
            "CD0": 0.0116,
        }

    def test_get_parameters(self):

        for key in self.params_baseline:
            param = self.ovl_solver.get_parameter(key)
            self.assertEqual(param, self.params_baseline[key], msg=key)

    def test_set_parameters(self):

        for key in self.params_baseline:
            # add each key to the update dict one at a time

            self.ovl_solver.set_parameter(key, self.params_baseline[key] + 0.1)

            # check that the parameter was updated
            param_updated = self.ovl_solver.get_parameter(key)
            self.assertEqual(param_updated, self.params_baseline[key] + 0.1, msg=key)

        # make sure the parameters stay set after an update
        self.ovl_solver.execute_run()
        for key in self.params_baseline:
            param_updated = self.ovl_solver.get_parameter(key)
            self.assertEqual(param_updated, self.params_baseline[key] + 0.1, msg=key)

    def test_set_fort_var(self):
        """
            test that the parameter changes effect the correct fortran variables
        """
        
        for key in self.params_baseline:
            # add each key to the update dict one at a time

            self.ovl_solver.set_parameter(key, self.params_baseline[key] + 0.1)

        # other parameters only get updated in exec subroutine
        self.ovl_solver.execute_run()
        
        self.assertEqual(self.ovl_solver.get_avl_fort_arr("CASE_R", "MACH"), 
                         self.params_baseline["Mach"] + 0.1)
        self.assertEqual(self.ovl_solver.get_avl_fort_arr("CASE_R", "CDREF"), 
                         self.params_baseline["CD0"] + 0.1)
        
        xyz_ref = self.ovl_solver.get_avl_fort_arr("CASE_R", "XYZREF")
        self.assertEqual(xyz_ref[0], self.params_baseline["X cg"] + 0.1)
        self.assertEqual(xyz_ref[1], self.params_baseline["Y cg"] + 0.1)
        self.assertEqual(xyz_ref[1], self.params_baseline["Z cg"] + 0.1)

class TestReferenceAPI(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.ref_data_baseline = {
            # "bank": 0.0,
            # "elevation": 0,
            "Sref": 1.5825, 
            "Cref": 0.35,
            "Bref": 4,
        }

    def test_get_ref_data(self):

        for key in self.ref_data_baseline:
            param = self.ovl_solver.get_reference_data()
            self.assertEqual(param[key], self.ref_data_baseline[key], msg=key)

    def test_set_ref_data(self):
        new_data = copy.deepcopy(self.ref_data_baseline)
        for key in self.ref_data_baseline:
            # add each key to the update dict one at a time
            new_data[key] += 0.1

            self.ovl_solver.set_reference_data(new_data)

            # check that the parameter was updated
            param_updated = self.ovl_solver.get_reference_data()
            self.assertEqual(param_updated[key], new_data[key], msg=key)

        # make sure the reference data stays set after an update
        self.ovl_solver.execute_run()
        for key in self.ref_data_baseline:
            param_updated = self.ovl_solver.get_reference_data()
            self.assertEqual(param_updated[key], new_data[key], msg=key)

if __name__ == "__main__":
    unittest.main()
    
 