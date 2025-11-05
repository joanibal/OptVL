"""This scipt demos looking at the airfoil data that is used by AVL"""

from optvl import OVLSolver
import matplotlib.pyplot as plt

ovl = OVLSolver(geo_file="../geom_files/aircraft.avl", debug=False)
surf_data = ovl.get_surface_params(include_geom=True, include_airfoils=True)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
idx_color = 0

for surf_key in ["Wing", "Horizontal Tail"]:
    # x coorindates for airfoil
    xlasec = surf_data[surf_key]["xlasec"]
    xuasec = surf_data[surf_key]["xuasec"]

    # upper and lower coordinates in the "thickness" direction
    zlasec = surf_data[surf_key]["zlasec"]
    zuasec = surf_data[surf_key]["zuasec"]

    # camber line
    casec = surf_data[surf_key]["casec"]

    # leading edge points
    yles = surf_data[surf_key]["yles"]
    xles = surf_data[surf_key]["xles"]

    # airfoil files
    afiles = surf_data[surf_key]["afiles"]

    for idx_airfoil in range(len(xlasec)):
        x_offset = xles[idx_airfoil]
        y_offset = yles[idx_airfoil]
        label = afiles[idx_airfoil]
        plt.plot(x_offset + xlasec[idx_airfoil], y_offset + zlasec[idx_airfoil], color=colors[idx_color], label=label)
        plt.plot(x_offset + xuasec[idx_airfoil], y_offset + zuasec[idx_airfoil], color=colors[idx_color])
        plt.plot(x_offset + xuasec[idx_airfoil], y_offset + casec[idx_airfoil], "--", color=colors[idx_color])
        idx_color += 1

plt.axis("equal")
plt.legend()
plt.show()
