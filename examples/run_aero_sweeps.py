from optvl import OVLSolver
import numpy as np

ovl = OVLSolver(geo_file="aircraft.avl", debug=False)

# set the angle of attack
ovl.set_constraint("alpha", 0.00)

# set the deflection of the elevator to trim the pitching moment
ovl.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")

ovl.set_parameter("Mach", 0.3)


print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    ovl.set_constraint("alpha", alpha)
    ovl.execute_run()
    run_data = ovl.get_total_forces()
    print(
        f" {alpha:10.6f}   {run_data['CL']:10.6f}   {run_data['CD']:10.6f}   {run_data['CDi']:10.6f}   {run_data['CDv']:10.6f}   {run_data['CM']:10.6f}"
    )

ovl.set_constraint("alpha", 0.00)

print("----------------- beta sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for beta in range(10):
    ovl.set_constraint("beta", beta)
    ovl.execute_run()
    run_data = ovl.get_total_forces()
    print(
        f" {beta:10.6f}   {run_data['CL']:10.6f}   {run_data['CD']:10.6f}   {run_data['CDi']:10.6f}   {run_data['CDv']:10.6f}   {run_data['CM']:10.6f}"
    )

ovl.set_constraint("beta", 0.00)

print("----------------- Mach sweep ----------------")
print("    Mach        Cl           Cd          Cdi          Cdv          Cm")
for mach in np.arange(0.0, 0.7, 0.1):
    ovl.set_parameter("Mach", mach)
    ovl.execute_run()
    run_data = ovl.get_total_forces()
    print(
        f" {mach:10.6f}   {run_data['CL']:10.6f}   {run_data['CD']:10.6f}   {run_data['CDi']:10.6f}   {run_data['CDv']:10.6f}   {run_data['CM']:10.6f}"
    )


print("----------------- CL sweep ----------------")
print("   Angle        Cl           Cd          Cdff          Cdv          Cm")
for cl in np.arange(0.6, 1.6, 0.1):
    ovl.set_trim_condition("CL", cl)
    ovl.execute_run()
    run_data = ovl.get_total_forces()
    alpha = ovl.get_parameter("alpha")
    print(
        f" {alpha:10.6f}   {run_data['CL']:10.6f}   {run_data['CD']:10.6f}   {run_data['CDi']:10.6f}   {run_data['CDv']:10.6f}   {run_data['CM']:10.6f}"
    )
