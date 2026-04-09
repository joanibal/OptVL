# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
import pickle


mesh = np.zeros((2,2,3))
mesh[1,:,0] = 1.0
mesh[:,1,1] = 1.0
mesh[:,:,2] = 1.0


surf = {
    "Wing": {
        # General
        # "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        # "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        # Geometry: Mesh
        "mesh": np.float64(mesh), # (nx,ny,3) numpy array containing mesh coordinates
        # Control Surface Specification
        "control_assignments": {
            "Elevator" : {"assignment":np.arange(0,mesh.shape[1]),
                      "xhinged": 0.5, # x/c location of hinge
                      "vhinged": np.array([0,1,0]), # vector giving hinge axis about which surface rotates
                      "gaind": -1.0, # control surface gain
                      "refld": 1.0  # control surface reflection, sign of deflection for duplicated surface
                      }
        },

    }
}

fuselage = {"Fuse pod": {
    # General
    # 'yduplicate': np.float64(0), # body is duplicated over the ysymm plane
    # Geometry
    "scale": np.array(
        [1.0, 1.0, 1.0]
    ),  # scaling factors applied to all x,y,z coordinates (chords areal so scaled by Xscale)
    "translate": np.array([0.0, 0.0, 0.0]),  # offset added on to all X,Y,Z values in this surface
    # Discretization
    "nvb": np.int32(4),  # number of source-line nodes
    "bspace": np.float64(2.0),  # lengthwise node spacing parameter
    "bfile": "../geom_files/fuseSimple.dat",  # body oml file name
}
}

input_dict = {
    "title": "MACH MDAO AVL",
    "mach": np.float64(0.0),
    "iysym": np.int32(0),
    "izsym": np.int32(0),
    "zsym": np.float64(0.0),
    "Sref": np.float64(1.123),
    "Cref": np.float64(0.25),
    "Bref": np.float64(6.01),
    "XYZref": np.array([0.0, 0, 0],dtype=np.float64),
    "CDp": np.float64(0.0),
    "surfaces": surf,
    "bodies": fuselage,
    # Global Control and DV info
    "dname": ["Elevator"],  # Name of control input for each corresonding index
}

# For verification
# base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
# geom_dir = os.path.join(base_dir, '..', 'geom_files')
# rect_file = os.path.join(geom_dir, 'rect_with_body.avl')


# solver = OVLSolver(input_dict=input_dict,debug=True)
# solver = OVLSolver(geo_file=rect_file,debug=True)

# solver.set_variable("alpha", 25.0)
# solver.set_variable("beta", 5.0)
# solver.execute_run()

# solver.plot_geom()

with open("rect_with_body.pkl", 'wb') as f:
    pickle.dump(input_dict, f)




