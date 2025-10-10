"""This scipt demos the use of the ways to vizualie the geometry and the solution in OptVL"""

from optvl import OVLSolver
import numpy as np
import matplotlib.pyplot as plt

ovl_solver = OVLSolver(geo_file="aircraft.avl", debug=False)
ovl_solver.plot_geom()

ovl_solver.set_constraint("alpha", 5.00)
ovl_solver.execute_run()

ovl_solver.plot_cp()
