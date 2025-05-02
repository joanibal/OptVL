from optvl import OVLSolver
import numpy as np

ovl_solver = OVLSolver(geo_file="aircraft.avl", debug=False)
alpha = 3.0
ovl_solver.set_constraint("alpha", alpha)
ovl_solver.execute_run()
run_data = ovl_solver.get_total_forces()
stab_derivs = ovl_solver.get_stab_derivs()
print('stability deriv', stab_derivs['dCL/dalpha'])

step = 1e-6
ovl_solver.set_constraint("alpha", alpha + step)
ovl_solver.execute_run()
run_data_p = ovl_solver.get_total_forces()
print("fd", (run_data_p["CL"] - run_data["CL"]) / step)
