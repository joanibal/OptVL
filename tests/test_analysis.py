# Extension modules
from optvl import OVLSolver
from ovl_avl_comparison_utils import run_comparison, set_inputs, get_avl_output_name, check_vals

# Standard Python modules
import os
import json

# External Python modules
import unittest
import numpy as np

# These test rely on reference data from avl
# create_avl_references.py is used to create this data
# modify or rerun it as needed


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "geom_files/supra.avl")
mass_file = os.path.join(base_dir, "geom_files/supra.mass")


class TestBasic(unittest.TestCase):
    """These tests are to check that a simple case works before launching the full matrix of tests used in the other tests"""

    def setUp(self):
        self.ovl = OVLSolver(geo_file=geom_file, mass_file=mass_file)

    def test_one_case(self):
        with open("avl_analysis_references/unconstrained_supra.json") as f:
            ref_data_cases = json.load(f)
            ref_data = ref_data_cases[0]

        # lets just look at the forces for one case
        set_inputs(self.ovl, ref_data)
        self.ovl.execute_run()
        force_data = self.ovl.get_total_forces()

        for key in force_data:
            avl_key = get_avl_output_name(key, self.ovl)
            avl_val = ref_data["outputs"]["total_forces"][avl_key]
            check_vals(force_data[key], avl_val, key, rtol=1e-15, printing=False)

class TestUnconstrained(unittest.TestCase):
    def test_aircraft(self):
        case = "aircraft"
        ovl = OVLSolver(
            geo_file=os.path.join(base_dir, f"geom_files/{case}.avl"),
            mass_file=os.path.join(base_dir, f"geom_files/{case}.mass"),
            timing=False,
            debug=False,
        )
        with open(f"avl_analysis_references/unconstrained_{case}.json") as f:
            ref_data_cases = json.load(f)
            run_comparison(ovl, ref_data_cases, printing=False)

    def test_supra(self):
        case = "supra"
        ovl = OVLSolver(
            geo_file=os.path.join(base_dir, f"geom_files/{case}.avl"),
            mass_file=os.path.join(base_dir, f"geom_files/{case}.mass"),
            timing=False,
            debug=False,
        )
        with open(f"avl_analysis_references/unconstrained_{case}.json") as f:
            ref_data_cases = json.load(f)
            run_comparison(ovl, ref_data_cases)


class TestConstrained(unittest.TestCase):
    def test_aircraft(self):
        case = "aircraft"
        ovl = OVLSolver(
            geo_file=os.path.join(base_dir, f"geom_files/{case}.avl"),
            mass_file=os.path.join(base_dir, f"geom_files/{case}.mass"),
            timing=False,
            debug=False,
        )
        with open(f"avl_analysis_references/unconstrained_{case}.json") as f:
            ref_data_cases = json.load(f)
            run_comparison(ovl, ref_data_cases)

    def test_supra(self):
        case = "supra"
        ovl = OVLSolver(
            geo_file=os.path.join(base_dir, f"geom_files/{case}.avl"),
            mass_file=os.path.join(base_dir, f"geom_files/{case}.mass"),
            timing=False,
            debug=False,
        )
        with open(f"avl_analysis_references/unconstrained_{case}.json") as f:
            ref_data_cases = json.load(f)
            run_comparison(ovl, ref_data_cases)


if __name__ == "__main__":
    unittest.main()
