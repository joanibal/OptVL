from optvl import OVLSolver
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt


ovl_solver = OVLSolver(geo_file="../geom_files/rectangle.avl", debug=False, timing=False)

# set the angle of attack
# ovl_solver.set_variable("alpha", 5.00)
yles = ovl_solver.get_surface_param("Wing", "yles")
xles = ovl_solver.get_surface_param("Wing", "xles")
zles = ovl_solver.get_surface_param("Wing", "zles")
span = yles[-1]
relative_span = yles / span


# Parameters
# params = ['xles', 'zles', 'yles']
values = {"xles": xles, "zles": zles, "yles": yles}

# Loop over each parameter type
for param, value in values.items():
    # Create subplots
    ax1 = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(2, 1, 2)
    ax2.set_ylabel("Z", rotation=0)
    ax2.set_xlabel("Y")
    ax1.set_ylabel("X", rotation=0)
    ax1.set_title(f"Modifying {param}")

    # Perform parameter sweep
    for d in np.linspace(-0.33 * span, 0.33 * span, 3):
        new_value = value + d * relative_span
        ovl_solver.set_surface_param("Wing", param, new_value)
        ovl_solver.plot_geom(axes=[ax1, ax2])

    # Reset to baseline
    ovl_solver.set_surface_param("Wing", param, value)

    # Show plot
    plt.show()
