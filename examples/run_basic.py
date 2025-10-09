from optvl import OVLSolver

ovl = OVLSolver(geo_file="aircraft.avl", debug=False)

# look at the geometry to see that everything is right
ovl.plot_geom()

# set the angle of attack
ovl.set_constraint("alpha", 0.00)

# set the deflection of the elevator to trim the pitching moment
ovl.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")

ovl.set_parameter("Mach", 0.3)

# This is the method that acutally runs the analysis
ovl.execute_run()

# print data about the run
force_data = ovl.get_total_forces()
print(f"CL:{force_data['CL']:10.6f}   CD:{force_data['CD']:10.6f}   CM:{force_data['CM']:10.6f}")

# lets look at the cp countours
ovl.plot_cp()
