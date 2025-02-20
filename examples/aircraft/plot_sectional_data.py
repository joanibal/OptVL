from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

ovl = AVLSolver(geo_file="aircraft.avl", debug=False)
ovl.set_constraint('alpha', 5.0)
ovl.execute_run()

strip_data = ovl.get_strip_forces()
for surf_key in strip_data:
    span_distance = strip_data[surf_key]['XYZ LE'][:,1]
    plt.plot(span_distance, strip_data[surf_key]['lift dist'], color='blue')
    plt.plot(span_distance, strip_data[surf_key]['CL'], color='red')
    plt.plot(span_distance, strip_data[surf_key]['CL perp'], color='firebrick', linestyle='--')

plt.legend(['lift dist', 'CL', 'CL perp.'])
plt.ylabel('lift distribution')
plt.xlabel('spanwise position')
plt.show()

