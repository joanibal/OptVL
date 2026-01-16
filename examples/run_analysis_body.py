from optvl import OVLSolver
import numpy as np

ovl = OVLSolver(geo_file="../geom_files/supra.avl", debug=False)
# set the angle of attack
ovl.set_parameter("Mach", 0.0)
ovl.execute_run()

body_names = ovl.get_body_names()
print("----------------- alpha sweep ----------------")
print("   Angle        CL           Cd        Cm        Cl")
for alpha in range(10):
    ovl.set_variable("alpha", alpha)
    ovl.execute_run()
    run_data = ovl.get_body_forces()
    for body_name in body_names:
        print(
            f' {alpha:10.2f}   {run_data[body_name]["CL"]:10.6e}   {run_data[body_name]["CD"]:10.6e}   {run_data[body_name]["Cm"]:10.6e}  {run_data[body_name]["Cl"]:10.6e}'
        )
