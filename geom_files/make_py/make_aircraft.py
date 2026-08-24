# =============================================================================
# Standard Python Modules
# =============================================================================
import os
import copy

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np
from openaerostruct.meshing.mesh_generator import generate_mesh
from openaerostruct.meshing.section_mesh_generator import generate_mesh as generate_section_mesh
from openaerostruct.geometry.utils import shear_z

# =============================================================================
# Local modules
# =============================================================================
from write_dict_to_py import write_dict_to_py
from optvl import OVLSolver


mesh_dict_wing = {
    "num_sections": 4,
    "sec_name": ["sec0", "sec1", "sec2", "sec3"],
    "symmetry": True,
    "ref_axis_pos": 0.0,
    "root_section": 0,
    # Geometry Parameters
    "taper": np.array([0.6666, 0.75, 0.8888, 1.0]),  # Wing taper for each section
    "span": np.array([0.25, 0.45, 0.55, 0.75]),  # Wing span for each section
    "sweep": np.array([0.0, 0.0, 0.0, 0,0]),  # Wing sweep for each section
    "root_chord": 0.45,  # Wing root chord for each section
    # Mesh Parameters
    "nx": 8,  # Number of chordwise panels. Same for all sections
    "ny": np.array([6, 6, 6, 6]),  # Number of spanwise panels for each section
}


_, sections_wing = generate_section_mesh(mesh_dict_wing)


mesh_dict_h_tail = {
        "num_y": 21,
        "num_x": 11,
        "wing_type": "rect",
        "symmetry": True,
        "span": 0.52*2,
        "root_chord": 0.325,
        "span_cos_spacing": 0.0,
        "chord_cos_spacing": 1.0,
    }
mesh_h_tail = generate_mesh(mesh_dict_h_tail)

mesh_dict_v_tail = {
        "num_y": 11,
        "num_x": 11,
        "wing_type": "rect",
        "symmetry": True,
        "span": 0.296*2,
        "root_chord": 0.325,
        "span_cos_spacing": 0.0,
        "chord_cos_spacing": 1.0,
    }
mesh_v_tail = generate_mesh(mesh_dict_v_tail)


shear_z(sections_wing[0],np.linspace(0.1,0.0,6))

for i in range(4):
    if i == 0:
        mesh_wing = sections_wing[i][:,:-1,:]
    elif i == 3:
        mesh_wing = np.hstack([mesh_wing,sections_wing[i]])
    else:
        mesh_wing = np.hstack([mesh_wing,sections_wing[i][:,:-1,:]])
# mesh_wing = np.hstack(sections_wing)

mesh_wing[:,:,1] = -mesh_wing[:,:,1]
mesh_wing = mesh_wing[:,::-1,:]

mesh_h_tail[:,:,1] = -mesh_h_tail[:,:,1]
mesh_h_tail = mesh_h_tail[:,::-1,:]

mesh_h_tail[:,:,0] += 1.60247523857

mesh_v_tail[:,:,2] = -mesh_v_tail[:,:,1]
mesh_v_tail[:,:,1] = 0.0

mesh_v_tail[:,:,0] += 1.60247523857


surf = {
    "Wing": {
        # General
        # "component": np.int32(0),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        # Geometry: Mesh
        "mesh": np.float64(mesh_wing), # (nx,ny,3) numpy array containing mesh coordinates
        # Airfoils
        # "afiles": np.array(["../geom_files/airfoils/A_1.dat","../geom_files/airfoils/A_2.dat","../geom_files/airfoils/A_3.dat","../geom_files/airfoils/A_4.dat","../geom_files/airfoils/A_5.dat"]),
        # "afiles": np.concatenate([np.repeat("../airfoils/A_1.dat",5),np.repeat("../airfoils/A_2.dat",5),np.repeat("../airfoils/A_3.dat",5),np.repeat("../airfoils/A_4.dat",5),["../airfoils/A_5.dat"]]),
        "afiles": "../geom_files/airfoils/ag40d.dat",

    },
    "Horizontal Tail": {
        # General
        # "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        "claf": 1.1078,
        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        # Geometry: Mesh
        "mesh": np.float64(mesh_h_tail), # (nx,ny,3) numpy array containing mesh coordinates
        # Airfoils
        "naca": "0012",
        # Control Surface Specification
        "control_assignments": {
            "Elevator" : {"assignment":np.arange(0,mesh_h_tail.shape[1]),
                      "xhinged": 0.0, # x/c location of hinge
                      "vhinged": np.array([0,1,0]), # vector giving hinge axis about which surface rotates
                      "gaind": -1.0, # control surface gain
                      "refld": 1.0  # control surface reflection, sign of deflection for duplicated surface
                      }
        },

    },
    "Vertical Tail": {
        # General
        "component": np.int32(2),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        # "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        "claf": 1.1078,
        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        # Geometry: Mesh
        "mesh": np.float64(mesh_v_tail), # (nx,ny,3) numpy array containing mesh coordinates
        # Airfoils
        "naca": "0012",
        # Control Surface Specification
        "control_assignments": {
            "Rudder" : {"assignment":np.arange(0,mesh_v_tail.shape[1]),
                      "xhinged": 0.5, # x/c location of hinge
                      "vhinged": np.array([0,0,1]), # vector giving hinge axis about which surface rotates
                      "gaind": 1.0, # control surface gain
                      "refld": -1.0  # control surface reflection, sign of deflection for duplicated surface
                      }
        },

    }
}


input_dict = {
    "title": "AIRCRAFT 1",
    "mach": np.float64(0.1),
    "iysym": np.int32(0),
    "izsym": np.int32(0),
    "zsym": np.float64(0.0),
    "Sref": np.float64(1.5825),
    "Cref": np.float64(0.35),
    "Bref": np.float64(4.0),
    "XYZref": np.array([0.084, 0.0, 0.0],dtype=np.float64),
    "CDp": np.float64(0.011600),
    "surfaces": surf,
    # "bodies": fuselage,
    # Global Control and DV info
    "dname": ["Elevator", "Rudder"],  # Name of control input for each corresonding index
}

# For verification
# base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
# geom_dir = os.path.join(base_dir, '..', 'geom_files')
# rect_file = os.path.join(geom_dir, 'aircraft.avl')


# solver = OVLSolver(input_dict=input_dict,debug=True)
# solver = OVLSolver(geo_file=rect_file,debug=True)

# solver.set_variable("alpha", 25.0)
# solver.set_variable("beta", 5.0)
# solver.execute_run()

# solver.plot_geom()
# solver.plot_cp()

# Write as a plain text Python file instead of pickle
output_file = "../aircraft.py"

write_dict_to_py(
    output_file,
    input_dict,
    var_name='input_dict',
    description='Generated by make_aircraft.py - Aircraft configuration with wing and tail'
)

print(f"Generated {output_file}")
