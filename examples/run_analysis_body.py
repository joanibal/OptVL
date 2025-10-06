from optvl import OVLSolver

ovl_solver = OVLSolver(geo_file="supra.avl", debug=False)
# set the angle of attack
ovl_solver.set_parameter("Mach", 0.0)

print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    ovl_solver.set_constraint("alpha", alpha)
    ovl_solver.execute_run()
    run_data = ovl_solver.get_total_forces()
    print(
        f" {alpha:10.6f}   {run_data['CL']:10.6f}   {run_data['CD']:10.6f}   {run_data['CDi']:10.6f}   {run_data['CDv']:10.6f}   {run_data['CM']:10.6f}"
    )
