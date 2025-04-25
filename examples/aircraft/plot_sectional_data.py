from optvl import AVLSolver
import numpy as np
import matplotlib.pyplot as plt

ovl = AVLSolver(geo_file="aircraft.avl", debug=False)
ovl.set_constraint('alpha', 5.0)
ovl.set_constraint('beta', 10.0)
ovl.execute_run()

# keys-start
strip_data = ovl.get_strip_forces()
first_surf = list(strip_data.keys())[0]
print(strip_data[first_surf].keys())
# keys-end

for surf_key in strip_data:
    span_distance = strip_data[surf_key]['Y LE']
    plt.plot(span_distance, strip_data[surf_key]['chord'], color='blue')
    plt.plot(span_distance, strip_data[surf_key]['twist'], color='red')

plt.legend(['chord', 'twist'])
plt.title('geometric spanwise data')
plt.xlabel('spanwise position')
plt.show()

for surf_key in strip_data:
    span_distance = strip_data[surf_key]['Y LE']
    plt.plot(span_distance, strip_data[surf_key]['lift dist'], color='blue')
    plt.plot(span_distance, strip_data[surf_key]['CL'], color='red')
    plt.plot(span_distance, strip_data[surf_key]['CL perp'], color='firebrick', linestyle='--')

plt.legend(['lift dist', 'CL', 'CL perp.'])
plt.title('lift spanwise data')
plt.xlabel('spanwise position')
plt.show()


strip_data = ovl.get_strip_forces()
for surf_key in strip_data:
    span_distance = strip_data[surf_key]['Y LE']
    plt.plot(span_distance, strip_data[surf_key]['CN'], color='C0')
    plt.plot(span_distance, strip_data[surf_key]['CR'], color='C1')

plt.legend(['roll distribution', 'yaw distribution'])
plt.title('roll and yaw spanwise data')
plt.xlabel('spanwise position')
plt.show()
