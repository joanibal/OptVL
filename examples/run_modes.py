from optvl import OVLSolver
import numpy as np
import matplotlib.pyplot as plt


def syssho(asys):
    """
    Prints out state-system matrix A in an organized manner.

    Parameters
    ----------
    asys : ndarray (nsys x nsys)
        State matrix A.
    """
    nsys = asys.shape[0]

    # Header
    state_labels = ["u", "w", "q", "the", "v", "p", "r", "phi", "x", "y", "z", "psi"]
    header = " ".join(f"{lab:>10}" for lab in state_labels) + "   |"
    print(header)

    # Rows
    for i in range(nsys):
        row_str = "".join(f"{val:11.4f}" for val in asys[i, :12])
        print(row_str)


ovl = OVLSolver(geo_file="aircraft.avl", mass_file="aircraft.mass", debug=False)

vel = 10.0

ovl.set_trim_condition("velocity", 10.0)
ovl.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")

ovl.execute_eigen_mode_calc()
Amat = ovl.get_system_matrix(in_body_axis=True)
syssho(Amat)


vals_ovl = ovl.get_eigenvalues()

# plot the eigenvalues
plt.plot(np.real(vals_ovl), np.imag(vals_ovl), "o")
plt.xlabel("real")
plt.ylabel("imag")
plt.title("Eigenvalues")
plt.grid("on")
plt.show()
# ----------------------------------------------------
