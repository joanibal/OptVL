from optvl import OVLSolver

ovl = OVLSolver(geo_file="../geom_files/aircraft.avl", debug=True)

# look at the geometry to see that everything is right
ovl.plot_geom()

# set the angle of attack
ovl.set_variable("alpha", 1.23)
# modify the mach number
ovl.set_parameter("Mach", 0.3)

# set the deflection of the elevator to trim the pitching moment
ovl.set_constraint("Elevator", "Cm", 0.00)

ovl.set_control_deflection("Elevator", 10.0)


# This is the method that acutally runs the analysis
ovl.execute_run()

# print data about the run 
force_data = ovl.get_total_forces()
print(
    f'CL:{force_data["CL"]:10.6f}   CD:{force_data["CD"]:10.6f}   Cm:{force_data["Cm"]:10.6f}'
)

# lets look at the cp countours 
ovl.plot_cp()