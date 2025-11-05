import openmdao.api as om
import numpy as np
import matplotlib.pyplot as plt
from optvl import OVLSolver
import glob
import os


cr = om.CaseReader("./run_opt_om_planform_out/opt_history.sql")
driver_cases = cr.list_cases("driver", out_stream=None)
obj_arr = np.zeros(len(driver_cases))
for idx_case in range(len(driver_cases)):
    obj_arr[idx_case] = cr.get_case(driver_cases[idx_case])["glide.duration"]

plt.plot(obj_arr)
plt.xlabel("iteration")
plt.ylabel("duration [s]")
plt.show()

# search the output directory for the latest file with a .dat extension
# Use glob to find all .dat files in the given directory
output_dir = "opt_output_sweep"
files = glob.glob(os.path.join(output_dir, "*.avl"))

# Find the .dat file with the latest modification time
latest_file = max(files, key=os.path.getmtime)

ovl = OVLSolver(geo_file=latest_file)

ovl.plot_geom()

ovl.execute_run()


strip_data = ovl.get_strip_forces()


# Create a figure and two subplots that share the x-axis
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))


for surf_key in strip_data:
    span_distance = strip_data[surf_key]["Y LE"]
    ax1.plot(span_distance, strip_data[surf_key]["chord"], color="blue", label="chord")
    ax1.plot(span_distance, strip_data[surf_key]["twist"], color="red", label="twist")
    ax2.plot(span_distance, strip_data[surf_key]["lift dist"], color="C0", label="list dist")

ax1.legend(["chord", "twist"])
ax2.legend(["lift dist"])
ax2.set_xlabel("span location")
plt.show()
