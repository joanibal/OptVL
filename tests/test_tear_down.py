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


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestCaseDerivs(unittest.TestCase):
    def setUp(self) -> None:
        self.ovl_solver = OVLSolver(geo_file=geom_file)

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_1(self):
        self.ovl_solver.execute_run()

    def test_2(self):
        self.ovl_solver.execute_run()


if __name__ == "__main__":
    unittest.main()
