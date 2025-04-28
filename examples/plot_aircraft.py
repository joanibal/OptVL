""" This scipt demos the use of the ways to vizualie the geometry and the solution in OptVL"""
from optvl import OVLSolver
import numpy as np
import matplotlib.pyplot as plt

avl_solver = OVLSolver(geo_file="aircraft.avl", debug=False)
avl_solver.plot_geom()

avl_solver.set_constraint("alpha", 5.00)
avl_solver.execute_run()

avl_solver.plot_cp()

