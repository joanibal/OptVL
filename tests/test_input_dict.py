# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

# =============================================================================
# Standard Python Modules
# =============================================================================
import os
from copy import deepcopy

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file1 = os.path.join(base_dir, "aircraft_wing_body.avl")
geom_file2 = os.path.join(base_dir, "aircraft_wing_body_naca.avl")
geom_file3 = os.path.join(base_dir, "aircraft_wing_body_airfoils.avl")
geom_file4 = os.path.join(base_dir, "aircraft_wing_body_afile.avl")


wing = {
        # General
        "num_sections": np.int32(2),
        "num_controls": np.array([1, 0], dtype=np.int32),
        "num_design_vars": np.array([1, 0], dtype=np.int32),
        "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        "wake": np.int32(
            1
        ),  # specifies that this surface is to NOT shed a wake, so that its strips will not have their Kutta conditions imposed
        "albe": np.int32(
            1
        ),  # specifies that this surface is unaffected by freestream direction changes specified by the alpha,beta angles and p,q,r rotation rates
        "load": np.int32(
            1
        ),  # specifies that the force and moment on this surface is to NOT be included in the overall forces and moments of the configuration
        "clcdsec": np.array(
            [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        ),  # profile-drag CD(CL) function for each section in this surface
        "cdcl": np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),  # profile-drag CD(CL) function for all sections in this surface, overrides Tahnks.
        "claf": np.array([1.0, 1.0]),  # CL alpha (dCL/da) scaling factor per section
        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        "angle": np.float64(0.0),  # offset added on to the Ainc values for all the defining sections in this surface
        "xles": np.array([0.0, 0.0]),  # leading edge cordinate vector(x component)
        "yles": np.array([-5.0, 0.0]),  # leading edge cordinate vector(y component)
        "zles": np.array([0.0, 0.0]),  # leading edge cordinate vector(z component)
        "chords": np.array([1.0, 1.0]),  # chord length vector
        "aincs": np.array([0.0, 0.0]),  # incidence angle vector
        # Geometry: Cross Sections
        "xfminmax": np.array([[0.0, 1.0], [0.0, 1.0]]),  # airfoil x/c limits
        # Paneling
        "nchordwise": np.int32(1),  # number of chordwise horseshoe vortice s placed on the surface
        "cspace": np.float64(0.0),  # chordwise vortex spacing parameter
        "nspan": np.int32(2),  # number of spanwise horseshoe vortices placed on the entire surface
        "sspace": np.float64(0.0),  # spanwise vortex spacing parameter for entire surface
        "nspans": np.array([5, 5], dtype=np.int32),  # number of spanwise elements vector
        "sspaces": np.array([3.0, 3.0], dtype=np.float64),  # spanwise spacing vector (for each section)
        "use surface spacing": np.int32(
            1
        ),  # surface spacing set under the surface heeading (known as LSURFSPACING in AVL)
    }

surf = {
    "Wing": {}
}

cont_surfs = {
     # Control Surfaces
    "icontd": [np.array([0],dtype=np.int32), np.array([],dtype=np.int32)],  # control variable index
    "xhinged": [np.array([0.8]), np.array([],dtype=np.float64)],  # x/c location of hinge
    "vhinged": [np.array([[0.0, 0.0, 0.0]]), np.array([],dtype=np.float64)],  # vector giving hinge axis about which surface rotates
    "gaind": [np.array([1.0]), np.array([],dtype=np.float64)],  # control surface gain
    "refld": [np.array([1.0]), np.array([],dtype=np.float64)],  # control surface reflection, sign of deflection for duplicated surface
}

des_var = {
    # Design Variables (AVL)
    "idestd": [np.array([0],dtype=np.int32), np.array([],dtype=np.int32)],  # design variable index
    "gaing": [np.array([1.0]), np.array([],dtype=np.float64)],  # desgin variable gain
}

cont_surf_names = {
    "dname": np.array(["flap"]),  # Name of control input for each corresonding index
}

des_var_names = {
    "gname": np.array(["des"]),  # Name of design var for each corresonding index
}

section_geom_naca = {
    # NACA
    'naca' : np.array(['2412','2412']), # 4-digit NACA airfoil
}

section_geom_raw = {
    # Coordinates
    'xasec': np.array([[0., 1.], [0., 1.]]), # the x coordinate aifoil section
    'casec': np.array([[0., 0.], [0., 0.]]), # camber line at xasec
    'tasec': np.array([[0., 0.], [0., 0.]]), # thickness at xasec
    'xuasec': np.array([[0., 0.], [0., 0.]]), # airfoil upper surface x-coords (alternative to specifying camber line)
    'xlasec': np.array([[0., 0.], [0., 0.]]),  # airfoil lower surface x-coords (alternative to specifying camber line)
    'zuasec': np.array([[0., 0.], [0., 0.]]),  # airfoil upper surface z-coords (alternative to specifying camber line)
    'zlasec': np.array([[0., 0.], [0., 0.]]),  # airfoil lower surface z-coords (alternative to specifying camber line)
}

section_geom_afiles = {
    # Airfoil Files
    'afiles': np.array([['ag40d.dat','ag40d.dat']]), # airfoil file names
}

section_geom_airfoils = {
    'airfoils': np.array([[[ 9.999910e-01,  9.941530e-01,  9.823030e-01,  9.688910e-01,
          9.552570e-01,  9.415770e-01,  9.279160e-01,  9.142470e-01,
          9.005640e-01,  8.868820e-01,  8.731920e-01,  8.595260e-01,
          8.458490e-01,  8.321950e-01,  8.185430e-01,  8.049090e-01,
          7.912780e-01,  7.776580e-01,  7.640420e-01,  7.588700e-01,
          7.542070e-01,  7.470050e-01,  7.364010e-01,  7.226770e-01,
          7.123140e-01,  6.950920e-01,  6.812120e-01,  6.672980e-01,
          6.560910e-01,  6.505180e-01,  6.393840e-01,  6.253940e-01,
          6.135510e-01,  5.974350e-01,  5.834900e-01,  5.695480e-01,
          5.584710e-01,  5.495080e-01,  5.414770e-01,  5.275530e-01,
          5.136300e-01,  4.996740e-01,  4.856880e-01,  4.716910e-01,
          4.576830e-01,  4.436490e-01,  4.295920e-01,  4.155110e-01,
          4.014190e-01,  3.873520e-01,  3.733220e-01,  3.593290e-01,
          3.453700e-01,  3.314470e-01,  3.175600e-01,  3.037050e-01,
          2.898840e-01,  2.760820e-01,  2.622900e-01,  2.485170e-01,
          2.347660e-01,  2.210300e-01,  2.073240e-01,  1.936520e-01,
          1.800080e-01,  1.664200e-01,  1.528790e-01,  1.394070e-01,
          1.260200e-01,  1.127260e-01,  9.956900e-02,  8.656800e-02,
          7.377700e-02,  6.123700e-02,  4.906700e-02,  3.745300e-02,
          2.677000e-02,  1.771000e-02,  1.097500e-02,  6.494000e-03,
          3.614000e-03,  1.761000e-03,  6.270000e-04,  8.200000e-05,
          1.700000e-05,  3.010000e-04,  1.155000e-03,  2.667000e-03,
          4.949000e-03,  8.338000e-03,  1.333800e-02,  2.024900e-02,
          2.909600e-02,  3.983400e-02,  5.177800e-02,  6.435200e-02,
          7.740100e-02,  9.104400e-02,  1.052410e-01,  1.197820e-01,
          1.344160e-01,  1.490330e-01,  1.636100e-01,  1.781480e-01,
          1.926230e-01,  2.070330e-01,  2.212630e-01,  2.354580e-01,
          2.496030e-01,  2.637340e-01,  2.778490e-01,  2.919720e-01,
          3.061040e-01,  3.202550e-01,  3.344330e-01,  3.486380e-01,
          3.628640e-01,  3.771310e-01,  3.914280e-01,  4.057110e-01,
          4.200090e-01,  4.342900e-01,  4.485450e-01,  4.627810e-01,
          4.769990e-01,  4.911820e-01,  5.053430e-01,  5.194980e-01,
          5.336440e-01,  5.477800e-01,  5.554190e-01,  5.651140e-01,
          5.755120e-01,  5.890550e-01,  6.025970e-01,  6.089950e-01,
          6.162710e-01,  6.299210e-01,  6.435540e-01,  6.533570e-01,
          6.623860e-01,  6.707840e-01,  6.844130e-01,  6.980280e-01,
          7.130020e-01,  7.248550e-01,  7.381830e-01,  7.438740e-01,
          7.479400e-01,  7.511780e-01,  7.535100e-01,  7.581350e-01,
          7.650660e-01,  7.785230e-01,  7.919490e-01,  8.053630e-01,
          8.189270e-01,  8.324880e-01,  8.460460e-01,  8.596250e-01,
          8.731990e-01,  8.867610e-01,  9.003580e-01,  9.139700e-01,
          9.276060e-01,  9.412600e-01,  9.549640e-01,  9.686680e-01,
          9.821570e-01,  9.940980e-01,  1.000009e+00],
        [ 4.780000e-04,  1.229000e-03,  2.910000e-03,  4.813000e-03,
          6.684000e-03,  8.567000e-03,  1.044600e-02,  1.236500e-02,
          1.429000e-02,  1.620600e-02,  1.810400e-02,  1.998200e-02,
          2.184500e-02,  2.369200e-02,  2.551500e-02,  2.731500e-02,
          2.908900e-02,  3.083400e-02,  3.254400e-02,  3.318400e-02,
          3.375700e-02,  3.463600e-02,  3.591900e-02,  3.753500e-02,
          3.872900e-02,  4.065000e-02,  4.213600e-02,  4.357500e-02,
          4.469800e-02,  4.524400e-02,  4.630800e-02,  4.760300e-02,
          4.866000e-02,  5.004300e-02,  5.118700e-02,  5.228200e-02,
          5.311800e-02,  5.377000e-02,  5.433700e-02,  5.527900e-02,
          5.616800e-02,  5.700700e-02,  5.779200e-02,  5.852000e-02,
          5.919000e-02,  5.980000e-02,  6.034900e-02,  6.083400e-02,
          6.125300e-02,  6.159900e-02,  6.187400e-02,  6.207300e-02,
          6.219200e-02,  6.222900e-02,  6.218100e-02,  6.204400e-02,
          6.181200e-02,  6.148200e-02,  6.104700e-02,  6.050200e-02,
          5.984000e-02,  5.905000e-02,  5.812500e-02,  5.706500e-02,
          5.585200e-02,  5.447500e-02,  5.292200e-02,  5.117600e-02,
          4.921700e-02,  4.702300e-02,  4.457400e-02,  4.183400e-02,
          3.878500e-02,  3.538400e-02,  3.160100e-02,  2.739000e-02,
          2.279400e-02,  1.807700e-02,  1.375300e-02,  1.017600e-02,
          7.287000e-03,  4.886000e-03,  2.821000e-03,  1.026000e-03,
         -4.880000e-04, -1.887000e-03, -3.319000e-03, -4.780000e-03,
         -6.332000e-03, -8.071000e-03, -1.002500e-02, -1.206600e-02,
         -1.403100e-02, -1.582400e-02, -1.730300e-02, -1.845100e-02,
         -1.932300e-02, -1.997500e-02, -2.043900e-02, -2.073900e-02,
         -2.090200e-02, -2.094200e-02, -2.088800e-02, -2.075600e-02,
         -2.055700e-02, -2.030500e-02, -2.001200e-02, -1.968400e-02,
         -1.932400e-02, -1.893500e-02, -1.851900e-02, -1.807800e-02,
         -1.761400e-02, -1.712700e-02, -1.662100e-02, -1.610000e-02,
         -1.556200e-02, -1.501300e-02, -1.445200e-02, -1.388200e-02,
         -1.330400e-02, -1.271900e-02, -1.213000e-02, -1.153500e-02,
         -1.093500e-02, -1.033300e-02, -9.728000e-03, -9.120000e-03,
         -8.512000e-03, -7.905000e-03, -7.576000e-03, -7.158000e-03,
         -6.712000e-03, -6.131000e-03, -5.553000e-03, -5.281000e-03,
         -4.972000e-03, -4.397000e-03, -3.826000e-03, -3.418000e-03,
         -3.042000e-03, -2.696000e-03, -2.138000e-03, -1.584000e-03,
         -9.810000e-04, -5.100000e-04,  1.600000e-05,  2.390000e-04,
          3.970000e-04,  4.810000e-04,  4.900000e-04,  5.060000e-04,
          5.310000e-04,  5.730000e-04,  6.050000e-04,  6.310000e-04,
          6.460000e-04,  6.510000e-04,  6.430000e-04,  6.220000e-04,
          5.880000e-04,  5.370000e-04,  4.690000e-04,  3.950000e-04,
          2.940000e-04,  1.800000e-04,  4.800000e-05, -9.000000e-05,
         -2.500000e-04, -4.030000e-04, -4.780000e-04]],

       [[ 9.999910e-01,  9.941530e-01,  9.823030e-01,  9.688910e-01,
          9.552570e-01,  9.415770e-01,  9.279160e-01,  9.142470e-01,
          9.005640e-01,  8.868820e-01,  8.731920e-01,  8.595260e-01,
          8.458490e-01,  8.321950e-01,  8.185430e-01,  8.049090e-01,
          7.912780e-01,  7.776580e-01,  7.640420e-01,  7.588700e-01,
          7.542070e-01,  7.470050e-01,  7.364010e-01,  7.226770e-01,
          7.123140e-01,  6.950920e-01,  6.812120e-01,  6.672980e-01,
          6.560910e-01,  6.505180e-01,  6.393840e-01,  6.253940e-01,
          6.135510e-01,  5.974350e-01,  5.834900e-01,  5.695480e-01,
          5.584710e-01,  5.495080e-01,  5.414770e-01,  5.275530e-01,
          5.136300e-01,  4.996740e-01,  4.856880e-01,  4.716910e-01,
          4.576830e-01,  4.436490e-01,  4.295920e-01,  4.155110e-01,
          4.014190e-01,  3.873520e-01,  3.733220e-01,  3.593290e-01,
          3.453700e-01,  3.314470e-01,  3.175600e-01,  3.037050e-01,
          2.898840e-01,  2.760820e-01,  2.622900e-01,  2.485170e-01,
          2.347660e-01,  2.210300e-01,  2.073240e-01,  1.936520e-01,
          1.800080e-01,  1.664200e-01,  1.528790e-01,  1.394070e-01,
          1.260200e-01,  1.127260e-01,  9.956900e-02,  8.656800e-02,
          7.377700e-02,  6.123700e-02,  4.906700e-02,  3.745300e-02,
          2.677000e-02,  1.771000e-02,  1.097500e-02,  6.494000e-03,
          3.614000e-03,  1.761000e-03,  6.270000e-04,  8.200000e-05,
          1.700000e-05,  3.010000e-04,  1.155000e-03,  2.667000e-03,
          4.949000e-03,  8.338000e-03,  1.333800e-02,  2.024900e-02,
          2.909600e-02,  3.983400e-02,  5.177800e-02,  6.435200e-02,
          7.740100e-02,  9.104400e-02,  1.052410e-01,  1.197820e-01,
          1.344160e-01,  1.490330e-01,  1.636100e-01,  1.781480e-01,
          1.926230e-01,  2.070330e-01,  2.212630e-01,  2.354580e-01,
          2.496030e-01,  2.637340e-01,  2.778490e-01,  2.919720e-01,
          3.061040e-01,  3.202550e-01,  3.344330e-01,  3.486380e-01,
          3.628640e-01,  3.771310e-01,  3.914280e-01,  4.057110e-01,
          4.200090e-01,  4.342900e-01,  4.485450e-01,  4.627810e-01,
          4.769990e-01,  4.911820e-01,  5.053430e-01,  5.194980e-01,
          5.336440e-01,  5.477800e-01,  5.554190e-01,  5.651140e-01,
          5.755120e-01,  5.890550e-01,  6.025970e-01,  6.089950e-01,
          6.162710e-01,  6.299210e-01,  6.435540e-01,  6.533570e-01,
          6.623860e-01,  6.707840e-01,  6.844130e-01,  6.980280e-01,
          7.130020e-01,  7.248550e-01,  7.381830e-01,  7.438740e-01,
          7.479400e-01,  7.511780e-01,  7.535100e-01,  7.581350e-01,
          7.650660e-01,  7.785230e-01,  7.919490e-01,  8.053630e-01,
          8.189270e-01,  8.324880e-01,  8.460460e-01,  8.596250e-01,
          8.731990e-01,  8.867610e-01,  9.003580e-01,  9.139700e-01,
          9.276060e-01,  9.412600e-01,  9.549640e-01,  9.686680e-01,
          9.821570e-01,  9.940980e-01,  1.000009e+00],
        [ 4.780000e-04,  1.229000e-03,  2.910000e-03,  4.813000e-03,
          6.684000e-03,  8.567000e-03,  1.044600e-02,  1.236500e-02,
          1.429000e-02,  1.620600e-02,  1.810400e-02,  1.998200e-02,
          2.184500e-02,  2.369200e-02,  2.551500e-02,  2.731500e-02,
          2.908900e-02,  3.083400e-02,  3.254400e-02,  3.318400e-02,
          3.375700e-02,  3.463600e-02,  3.591900e-02,  3.753500e-02,
          3.872900e-02,  4.065000e-02,  4.213600e-02,  4.357500e-02,
          4.469800e-02,  4.524400e-02,  4.630800e-02,  4.760300e-02,
          4.866000e-02,  5.004300e-02,  5.118700e-02,  5.228200e-02,
          5.311800e-02,  5.377000e-02,  5.433700e-02,  5.527900e-02,
          5.616800e-02,  5.700700e-02,  5.779200e-02,  5.852000e-02,
          5.919000e-02,  5.980000e-02,  6.034900e-02,  6.083400e-02,
          6.125300e-02,  6.159900e-02,  6.187400e-02,  6.207300e-02,
          6.219200e-02,  6.222900e-02,  6.218100e-02,  6.204400e-02,
          6.181200e-02,  6.148200e-02,  6.104700e-02,  6.050200e-02,
          5.984000e-02,  5.905000e-02,  5.812500e-02,  5.706500e-02,
          5.585200e-02,  5.447500e-02,  5.292200e-02,  5.117600e-02,
          4.921700e-02,  4.702300e-02,  4.457400e-02,  4.183400e-02,
          3.878500e-02,  3.538400e-02,  3.160100e-02,  2.739000e-02,
          2.279400e-02,  1.807700e-02,  1.375300e-02,  1.017600e-02,
          7.287000e-03,  4.886000e-03,  2.821000e-03,  1.026000e-03,
         -4.880000e-04, -1.887000e-03, -3.319000e-03, -4.780000e-03,
         -6.332000e-03, -8.071000e-03, -1.002500e-02, -1.206600e-02,
         -1.403100e-02, -1.582400e-02, -1.730300e-02, -1.845100e-02,
         -1.932300e-02, -1.997500e-02, -2.043900e-02, -2.073900e-02,
         -2.090200e-02, -2.094200e-02, -2.088800e-02, -2.075600e-02,
         -2.055700e-02, -2.030500e-02, -2.001200e-02, -1.968400e-02,
         -1.932400e-02, -1.893500e-02, -1.851900e-02, -1.807800e-02,
         -1.761400e-02, -1.712700e-02, -1.662100e-02, -1.610000e-02,
         -1.556200e-02, -1.501300e-02, -1.445200e-02, -1.388200e-02,
         -1.330400e-02, -1.271900e-02, -1.213000e-02, -1.153500e-02,
         -1.093500e-02, -1.033300e-02, -9.728000e-03, -9.120000e-03,
         -8.512000e-03, -7.905000e-03, -7.576000e-03, -7.158000e-03,
         -6.712000e-03, -6.131000e-03, -5.553000e-03, -5.281000e-03,
         -4.972000e-03, -4.397000e-03, -3.826000e-03, -3.418000e-03,
         -3.042000e-03, -2.696000e-03, -2.138000e-03, -1.584000e-03,
         -9.810000e-04, -5.100000e-04,  1.600000e-05,  2.390000e-04,
          3.970000e-04,  4.810000e-04,  4.900000e-04,  5.060000e-04,
          5.310000e-04,  5.730000e-04,  6.050000e-04,  6.310000e-04,
          6.460000e-04,  6.510000e-04,  6.430000e-04,  6.220000e-04,
          5.880000e-04,  5.370000e-04,  4.690000e-04,  3.950000e-04,
          2.940000e-04,  1.800000e-04,  4.800000e-05, -9.000000e-05,
         -2.500000e-04, -4.030000e-04, -4.780000e-04]]]),
}

fuselage = {
        # General
        # 'yduplicate': np.float64(0), # body is duplicated over the ysymm plane
        # Geometry
        "scale": np.array(
            [[1.0, 1.0, 1.0]]
        ),  # scaling factors applied to all x,y,z coordinates (chords areal so scaled by Xscale)
        "translate": np.array([[0.0, 0.0, 0.0]]),  # offset added on to all X,Y,Z values in this surface
        # Discretization
        "nvb": np.int32(10),  # number of source-line nodes
        "bspace": np.float64(0.0),  # lengthwise node spacing parameter
    }

body = {"Fuselage": {}}

body_oml_file = {
    "bfile": "fuseSupra.dat",  # body oml file name
}

body_oml = {
    'body_oml': np.array([[ 5.1000000e+01,  1.3000000e+01,  1.2501041e+01,  1.1990128e+01,
         1.1478845e+01,  1.0967545e+01,  1.0456449e+01,  9.9455080e+00,
         9.4346670e+00,  8.9242760e+00,  8.4146700e+00,  7.9052190e+00,
         7.3938250e+00,  6.8815310e+00,  6.3708380e+00,  5.8623540e+00,
         5.3552870e+00,  4.8516350e+00,  4.3529120e+00,  3.8587060e+00,
         3.3684150e+00,  2.8818080e+00,  2.3989040e+00,  1.9199830e+00,
         1.4443960e+00,  9.7026900e-01,  4.9432400e-01,  1.3164000e-02,
        -4.6862800e-01, -9.3586500e-01, -1.3963740e+00, -1.8608640e+00,
        -2.3308510e+00, -2.8055690e+00, -3.2844620e+00, -3.7678000e+00,
        -4.2559290e+00, -4.7485880e+00, -5.2449740e+00, -5.7440110e+00,
        -6.2444730e+00, -6.7458890e+00, -7.2478060e+00, -7.7498200e+00,
        -8.2515490e+00, -8.7527420e+00, -9.2531120e+00, -9.7524140e+00,
        -1.0250331e+01, -1.0746588e+01, -1.1240778e+01, -1.1732254e+01,
        -1.2220179e+01, -1.2703416e+01, -1.3180237e+01, -1.3647286e+01,
        -1.4098492e+01, -1.4522828e+01, -1.4903724e+01, -1.5225227e+01,
        -1.5482801e+01, -1.5684700e+01, -1.5843814e+01, -1.5970662e+01,
        -1.6074086e+01, -1.6160914e+01, -1.6234984e+01, -1.6298600e+01,
        -1.6353098e+01, -1.6398832e+01, -1.6436125e+01, -1.6464678e+01,
        -1.6484758e+01, -1.6496397e+01, -1.6500000e+01, -1.6496397e+01,
        -1.6484758e+01, -1.6464678e+01, -1.6436125e+01, -1.6398832e+01,
        -1.6353098e+01, -1.6298600e+01, -1.6234984e+01, -1.6160914e+01,
        -1.6074086e+01, -1.5970662e+01, -1.5843814e+01, -1.5684700e+01,
        -1.5482801e+01, -1.5225227e+01, -1.4903724e+01, -1.4522828e+01,
        -1.4098492e+01, -1.3647286e+01, -1.3180237e+01, -1.2703416e+01,
        -1.2220179e+01, -1.1732254e+01, -1.1240778e+01, -1.0746588e+01,
        -1.0250331e+01, -9.7524140e+00, -9.2531120e+00, -8.7527420e+00,
        -8.2515490e+00, -7.7498200e+00, -7.2478060e+00, -6.7458890e+00,
        -6.2444730e+00, -5.7440110e+00, -5.2449740e+00, -4.7485880e+00,
        -4.2559290e+00, -3.7678000e+00, -3.2844620e+00, -2.8055690e+00,
        -2.3308510e+00, -1.8608640e+00, -1.3963740e+00, -9.3586500e-01,
        -4.6862800e-01,  1.3164000e-02,  4.9432400e-01,  9.7026900e-01,
         1.4443960e+00,  1.9199830e+00,  2.3989040e+00,  2.8818080e+00,
         3.3684150e+00,  3.8587060e+00,  4.3529120e+00,  4.8516350e+00,
         5.3552870e+00,  5.8623540e+00,  6.3708380e+00,  6.8815310e+00,
         7.3938250e+00,  7.9052190e+00,  8.4146700e+00,  8.9242760e+00,
         9.4346670e+00,  9.9455080e+00,  1.0456449e+01,  1.0967545e+01,
         1.1478845e+01,  1.1990128e+01,  1.2501041e+01,  1.3000000e+01,
         1.3000000e+01,  5.1000000e+01],
       [ 2.5000000e-01,  4.2500000e-01,  4.3400200e-01,  4.4378100e-01,
         4.5416300e-01,  4.6516600e-01,  4.7681600e-01,  4.8916400e-01,
         5.0223000e-01,  5.1601600e-01,  5.3052100e-01,  5.4577300e-01,
         5.6185000e-01,  5.7873200e-01,  5.9634800e-01,  6.1468700e-01,
         6.3380400e-01,  6.5366500e-01,  6.7427100e-01,  6.9575800e-01,
         7.1831600e-01,  7.4211300e-01,  7.6726300e-01,  7.9379300e-01,
         8.2157000e-01,  8.5053800e-01,  8.8048300e-01,  9.1093000e-01,
         9.4074300e-01,  9.6817700e-01,  9.9301600e-01,  1.0152710e+00,
         1.0347150e+00,  1.0512260e+00,  1.0648550e+00,  1.0757730e+00,
         1.0841740e+00,  1.0901930e+00,  1.0939120e+00,  1.0953920e+00,
         1.0946750e+00,  1.0917780e+00,  1.0866900e+00,  1.0793810e+00,
         1.0698160e+00,  1.0579560e+00,  1.0437520e+00,  1.0271300e+00,
         1.0080050e+00,  9.8621100e-01,  9.6157400e-01,  9.3389600e-01,
         9.0294400e-01,  8.6843000e-01,  8.3000600e-01,  7.8738500e-01,
         7.4051400e-01,  6.9001900e-01,  6.3778800e-01,  5.8680700e-01,
         5.3960100e-01,  4.9702400e-01,  4.5858000e-01,  4.2358200e-01,
         3.9094200e-01,  3.5940800e-01,  3.2802400e-01,  2.9559300e-01,
         2.6097800e-01,  2.2352300e-01,  1.8313000e-01,  1.3933200e-01,
         9.3422000e-02,  4.6702000e-02,  0.0000000e+00, -4.6702000e-02,
        -9.3422000e-02, -1.3933200e-01, -1.8313000e-01, -2.2352300e-01,
        -2.6097800e-01, -2.9559300e-01, -3.2802400e-01, -3.5940800e-01,
        -3.9094200e-01, -4.2358200e-01, -4.5858000e-01, -4.9702400e-01,
        -5.3960100e-01, -5.8680700e-01, -6.3778800e-01, -6.9001900e-01,
        -7.4051400e-01, -7.8738500e-01, -8.3000600e-01, -8.6843000e-01,
        -9.0294400e-01, -9.3389600e-01, -9.6157400e-01, -9.8621100e-01,
        -1.0080050e+00, -1.0271300e+00, -1.0437520e+00, -1.0579560e+00,
        -1.0698160e+00, -1.0793810e+00, -1.0866900e+00, -1.0917780e+00,
        -1.0946750e+00, -1.0953920e+00, -1.0939120e+00, -1.0901930e+00,
        -1.0841740e+00, -1.0757730e+00, -1.0648550e+00, -1.0512260e+00,
        -1.0347150e+00, -1.0152710e+00, -9.9301600e-01, -9.6817700e-01,
        -9.4074300e-01, -9.1093000e-01, -8.8048300e-01, -8.5053800e-01,
        -8.2157000e-01, -7.9379300e-01, -7.6726300e-01, -7.4211300e-01,
        -7.1831600e-01, -6.9575800e-01, -6.7427100e-01, -6.5366500e-01,
        -6.3380400e-01, -6.1468700e-01, -5.9634800e-01, -5.7873200e-01,
        -5.6185000e-01, -5.4577300e-01, -5.3052100e-01, -5.1601600e-01,
        -5.0223000e-01, -4.8916400e-01, -4.7681600e-01, -4.6516600e-01,
        -4.5416300e-01, -4.4378100e-01, -4.3400200e-01, -4.2500000e-01,
        -4.2500000e-01, -2.5000000e-01]]) # raw body oml coords
}


input_dict = {
    "title": "Aircraft",
    "mach": np.float64(0.0),
    "iysym": np.int32(0),
    "izsym": np.int32(0),
    "zsymm": np.float64(0.0),
    "Sref": np.float64(10.0),
    "Cref": np.float64(1.0),
    "Bref": np.float64(10.0),
    "Xref": np.float64(0.25),
    "Yref": np.float64(0.0),
    "Zref": np.float64(0.0),
    "CDp": np.float64(0.0),
    "surfaces": {},
    "bodies": {},
}




class TestGeom(unittest.TestCase):
    def setUp(self):
        # Setup all 5 cases for testing

        # Flat Camberline
        input_dict1 = deepcopy(input_dict)
        surf1 = deepcopy(surf)
        surf1["Wing"] = wing|cont_surfs|des_var
        body1 = deepcopy(body)
        body1["Fuselage"] = fuselage|body_oml
        input_dict1 = input_dict1|cont_surf_names|des_var_names
        input_dict1["surfaces"] = surf1
        input_dict1["bodies"] = body1

        # NACA 2412
        input_dict2 = deepcopy(input_dict)
        surf2 = deepcopy(surf)
        surf2["Wing"] = wing|section_geom_naca|cont_surfs|des_var
        body1 = deepcopy(body)
        body1["Fuselage"] = fuselage|body_oml
        input_dict2 = input_dict2|cont_surf_names|des_var_names
        input_dict2["surfaces"] = surf2
        input_dict2["bodies"] = body1

        # Airfoil Coordinates
        input_dict3 = deepcopy(input_dict)
        surf3 = deepcopy(surf)
        surf3["Wing"] = wing|section_geom_airfoils|cont_surfs|des_var
        body1 = deepcopy(body)
        body1["Fuselage"] = fuselage|body_oml
        input_dict3 = input_dict3|cont_surf_names|des_var_names
        input_dict3["surfaces"] = surf3
        input_dict3["bodies"] = body1

        # Airfoil Files
        input_dict4 = deepcopy(input_dict)
        surf4 = deepcopy(surf)
        surf4["Wing"] = wing|section_geom_airfoils|cont_surfs|des_var
        body1 = deepcopy(body)
        body1["Fuselage"] = fuselage|body_oml
        input_dict4 = input_dict4|cont_surf_names|des_var_names
        input_dict4["surfaces"] = surf4
        input_dict4["bodies"] = body1

        # Body Files + Airfoil Files
        input_dict5 = deepcopy(input_dict)
        surf5 = deepcopy(surf)
        surf5["Wing"] = wing|section_geom_airfoils|cont_surfs|des_var
        body2 = deepcopy(body)
        body2["Fuselage"] = fuselage|body_oml_file
        input_dict5 = input_dict5|cont_surf_names|des_var_names
        input_dict5["surfaces"] = surf5
        input_dict5["bodies"] = body2


        # Solvers loaded with inputs dicts
        self.ovl_solver_1 = OVLSolver(input_dict=input_dict1)
        self.ovl_solver_2 = OVLSolver(input_dict=input_dict2)
        self.ovl_solver_3 = OVLSolver(input_dict=input_dict3)
        self.ovl_solver_4 = OVLSolver(input_dict=input_dict4)
        self.ovl_solver_5 = OVLSolver(input_dict=input_dict5)

        #TODO : Fix file writing
        # Solvers loaded with inputs files (case 4 and 5 use the same file)
        self.ovl_solver_f1 = OVLSolver(geo_file=geom_file1)
        self.ovl_solver_f2 = OVLSolver(geo_file=geom_file2)
        self.ovl_solver_f3 = OVLSolver(geo_file=geom_file3)
        self.ovl_solver_f4 = OVLSolver(geo_file=geom_file4)

        self.solvers = [self.ovl_solver_1,self.ovl_solver_2,self.ovl_solver_3,self.ovl_solver_4,self.ovl_solver_5]
        self.solvers_f = [self.ovl_solver_f1,self.ovl_solver_f2,self.ovl_solver_f3,self.ovl_solver_f4]

    def test_surface_params(self):
        data = []

        for solver in self.solvers:
            data.append(solver.get_surface_params(include_geom=True, include_paneling=True, include_con_surf=True))

        from pprint import pprint

        for i in range(len(self.solvers)):
            for s in surf:
                for key in surf[s]:
                    np.testing.assert_allclose(
                        data[i][s][key],
                        surf[s][key],
                        atol=1e-8,
                        err_msg=f"Surface `{s}` key `{key}` does not match reference data",
                    )

        # self.ovl_solver.set_constraint("alpha", 6.00)
        # self.ovl_solver.set_constraint("beta", 2.00)
        # self.ovl_solver.execute_run()
        
        # assert self.ovl_solver.get_num_surfaces() == 5
        # assert self.ovl_solver.get_num_strips() == 90
        # assert self.ovl_solver.get_mesh_size() == 780

        # np.testing.assert_allclose(
        #     self.ovl_solver.get_constraint("alpha"),
        #     6.0,
        #     rtol=1e-8,
        # )
        # np.testing.assert_allclose(
        #     self.ovl_solver.get_constraint("beta"),
        #     2.0,
        #     rtol=1e-8,
        # )
        
        # coefs = self.ovl_solver.get_total_forces()
        # np.testing.assert_allclose(
        #     coefs["CL"],
        #     5.407351081559913,
        #     rtol=1e-8,
        # )

        # self.ovl_solver.set_surface_params(data)
        

        # assert self.ovl_solver.get_num_surfaces() == 5
        # assert self.ovl_solver.get_num_strips() == 90
        # assert self.ovl_solver.get_mesh_size() == 780

        # self.ovl_solver.set_constraint("alpha", 6.00)
        # self.ovl_solver.set_constraint("beta", 2.00)
        # self.ovl_solver.execute_run()

        # np.testing.assert_allclose(
        #     self.ovl_solver.get_constraint("alpha"),
        #     6.0,
        #     rtol=1e-8,
        # )
        # np.testing.assert_allclose(
        #     self.ovl_solver.get_constraint("beta"),
        #     2.0,
        #     rtol=1e-8,
        # )
        
        # coefs = self.ovl_solver.get_total_forces()
        # np.testing.assert_allclose(
        #     coefs["CL"],
        #     5.407351081559913,
        #     rtol=1e-8,
        # )



if __name__ == "__main__":
    unittest.main()
