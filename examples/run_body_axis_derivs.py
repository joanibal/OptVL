from optvl import OVLSolver

ovl_solver = OVLSolver(geo_file="../geom_files/aircraft.avl", debug=False)
alpha = 3.0
ovl_solver.set_variable("alpha", alpha)
ovl_solver.execute_run()

derivs = ovl_solver.get_body_axis_derivs()

for key in derivs:
    print(f"{key}: {derivs[key]}")