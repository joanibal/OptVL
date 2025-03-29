from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)
avl_solver.set_constraint("alpha", 5.00)
avl_solver.execute_run()
avl_solver.write_tecplot('test_file')