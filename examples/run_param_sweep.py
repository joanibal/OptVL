from optvl import OVLSolver
import numpy as np

write_tecplot_files = True

ovl_solver = OVLSolver(geo_file="aircraft.avl", debug=False, timing=False)

# set the angle of attack
ovl_solver.set_constraint("alpha", 5.00)

for idx_scale, y_scale in enumerate(np.linspace(0.5, 1.5, 5)):
    ovl_solver.set_surface_params({"Wing": {"scale": np.array([1, y_scale, 1])}})

    ovl_solver.execute_run()
    stab_derivs = ovl_solver.get_stab_derivs()

    print(f"----------------- y_scale: {y_scale} ----------------")
    for key in stab_derivs:
        print(f"{key:16}: {stab_derivs[key]:.6f}")

    if write_tecplot_files:
        # this way works on tecplot and paraview
        ovl_solver.write_tecplot(f"wing_scale_{idx_scale}")

        # Warning: The solution time does not work on paraview
        # ovl_solver.write_tecplot(f'wing_scale_{y_scale}', solution_time=idx_scale)
