from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

avl = AVLSolver(geo_file="aircraft.avl", mass_file="aircraft.mass",  debug=False)

vel = 10
avl.set_parameter("velocity", vel)
dens = avl.get_parameter("density")
g = avl.get_parameter("grav.acc.")
mass = avl.get_parameter("mass")
weight = mass * g
cl = weight / (0.5 * dens * vel**2)
avl.set_trim_condition("CL", cl)

avl.execute_eigen_mode_calc()

vals_avl = avl.get_eigenvalues()

# plot the eigenvalues
plt.plot(np.real(vals_avl),np.imag(vals_avl), 'o')
plt.xlabel('real')
plt.ylabel('imag')
plt.title('Eigenvalues')
plt.show()
# ----------------------------------------------------
