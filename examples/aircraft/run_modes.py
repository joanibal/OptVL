from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

ovl = AVLSolver(geo_file="aircraft.avl", mass_file="aircraft.mass",  debug=False)

vel = 10
ovl.set_parameter("velocity", vel)
dens = ovl.get_parameter("density")
g = ovl.get_parameter("grav.acc.")
mass = ovl.get_parameter("mass")
weight = mass * g
cl = weight / (0.5 * dens * vel**2)
ovl.set_trim_condition("CL", cl)
ovl.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")

ovl.execute_eigen_mode_calc()

vals_ovl = ovl.get_eigenvalues()

# plot the eigenvalues
plt.plot(np.real(vals_ovl),np.imag(vals_ovl), 'o')
plt.xlabel('real')
plt.ylabel('imag')
plt.title('Eigenvalues')
plt.grid('on')
plt.show()
# ----------------------------------------------------
