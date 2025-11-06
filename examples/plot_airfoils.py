"""This scipt demos looking at the airfoil data that is used by AVL"""

from optvl import OVLSolver
import matplotlib.pyplot as plt

ovl = OVLSolver(geo_file="../geom_files/aircraft.avl", debug=True)
surf_data = ovl.get_surface_params(include_geom=True, include_airfoils=True)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
idx_color = 0

for surf_key in ["Wing", "Horizontal Tail"]:
    # x coorindates for airfoil
    airfoil_coords = surf_data[surf_key]["airfoils"]
    
    # leading edge points
    yles = surf_data[surf_key]["yles"]
    xles = surf_data[surf_key]["xles"]

    # airfoil files
    afiles = surf_data[surf_key]["afiles"]

    for idx_airfoil in range(len(airfoil_coords)):
        x_offset = xles[idx_airfoil]
        y_offset = yles[idx_airfoil]
        coords = airfoil_coords[idx_airfoil]
        
        label = afiles[idx_airfoil]
        plt.plot(x_offset + coords[0,:], y_offset + coords[1,:], color=colors[idx_color], label=label)
        idx_color += 1

plt.axis("equal")
plt.legend()
plt.show()
