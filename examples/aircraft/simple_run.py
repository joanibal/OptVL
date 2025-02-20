from optvl import AVLSolver
import numpy as np

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)

# set the angle of attack
avl_solver.set_constraint("alpha", 0.00)

# control surface names from geometry file
avl_solver.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")
avl_solver.set_constraint("Rudder", 0.00, con_var="Cn yaw moment")

avl_solver.set_parameter("Mach", 0.3)

# This is the method that acutally runs the analysis
avl_solver.execute_run()

print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    avl_solver.set_constraint("alpha", alpha)
    avl_solver.execute_run()
    run_data = avl_solver.get_total_forces()
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )


print("----------------- CL sweep ----------------")
print("   Angle        Cl           Cd           Cdff         Cdv          Cm         CN")
for cl in np.arange(0.6, 1.7, 0.1):
    avl_solver.set_trim_condition("CL", cl)
    avl_solver.execute_run()
    run_data = avl_solver.get_total_forces()
    alpha = avl_solver.get_parameter("alpha")
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f} {run_data["CN SA"]:10.6f} '
    )
