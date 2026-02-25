# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np



mesh = np.load("wing_mesh.npy")

surf = {
    "Wing": {
        # General
        "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        # "wake": np.int32(
        #     1
        # ),  # specifies that this surface is to NOT shed a wake, so that its strips will not have their Kutta conditions imposed
        # "albe": np.int32(
        #     1
        # ),  # specifies that this surface is unaffected by freestream direction changes specified by the alpha,beta angles and p,q,r rotation rates
        # "load": np.int32(
        #     1
        # ),  # specifies that the force and moment on this surface is to NOT be included in the overall forces and moments of the configuration
        # "clcdsec": np.array(
        #     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        # ),  # profile-drag CD(CL) function for each section in this surface (provide a single entry and OptVL applies to all strips, otherwise provide a vector corresponding to each strip)
        # "cdcl": np.array(
        #     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        # ),  # profile-drag CD(CL) function for all sections in this surface, overrides Tahnks.
        "claf": 1.0,  # CL alpha (dCL/da) scaling factor per section (provide a single entry and OptVL applies to all strips, otherwise provide a vector corresponding to each strip)

        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        "angle": np.float64(0.0),  # offset added on to the Ainc values for all the defining sections in this surface
        "aincs": np.ones(mesh.shape[1]), # incidence angle vector (provide a single entry and OptVL applies to all strips, otherwise provide a vector corresponding to each strip)

        # Geometry: Mesh
        "mesh": np.float64(mesh), # (nx,ny,3) numpy array containing mesh coordinates
        "flatten mesh": True, # True by default so can be turned off or just excluded (not recommended)

        # Geometry: Cross Sections (provide a single entry and OptVL applies to all strips, otherwise provide a vector corresponding to each strip)
        # "xfminmax": np.array([[0.0, 1.0]]),  # airfoil x/c limits
        # NACA
        # 'naca' : '2412', # 4-digit NACA airfoil
        # Direct Assignment of camberline/thickness
        # 'xasec': np.array([[0., 1.]]), # the x coordinate aifoil section
        # 'casec': np.array([[0., 0.]]), # camber line at xasec
        # 'tasec': np.array([[0., 0.]]), # thickness at xasec
        # 'xuasec': np.array([[0., 0.]]), # airfoil upper surface x-coords (alternative to specifying camber line)
        # 'xlasec': np.array([[0., 0.]]),  # airfoil lower surface x-coords (alternative to specifying camber line)
        # 'zuasec': np.array([[0., 0.]]),  # airfoil upper surface z-coords (alternative to specifying camber line)
        # 'zlasec': np.array([[0., 0.]]),  # airfoil lower surface z-coords (alternative to specifying camber line)
        # Airfoil Files
        'afiles': 'airfoils/ag40d.dat', # airfoil file names


        # Control Surface Specification
        "control_assignments": {
            "flap" : {"assignment":np.arange(0,mesh.shape[1]),
                      "xhinged": 0.8, # x/c location of hinge
                      "vhinged": np.zeros(3), # vector giving hinge axis about which surface rotates
                      "gaind": 1.0, # control surface gain
                      "refld": 1.0  # control surface reflection, sign of deflection for duplicated surface
                      }
        },

        # Design Variables (AVL) Specification
        "design_var_assignments": {
            "des" : {"assignment":np.arange(0,mesh.shape[1]),
                     "gaing":1.0}
        },
    }
}


surf_avl = {
    "Wing": {
        # General
        "num_sections": np.int32(2),
        "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        # "wake": np.int32(
        #     1
        # ),  # specifies that this surface is to NOT shed a wake, so that its strips will not have their Kutta conditions imposed
        # "albe": np.int32(
        #     1
        # ),  # specifies that this surface is unaffected by freestream direction changes specified by the alpha,beta angles and p,q,r rotation rates
        # "load": np.int32(
        #     1
        # ),  # specifies that the force and moment on this surface is to NOT be included in the overall forces and moments of the configuration
        # "clcdsec": np.array(
        #     [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        # ),  # profile-drag CD(CL) function for each section in this surface
        # "cdcl": np.array(
        #     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        # ),  # profile-drag CD(CL) function for all sections in this surface, overrides Tahnks.
        "claf": np.array([1.0, 1.0]),  # CL alpha (dCL/da) scaling factor per section

        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        # "angle": np.float64(0.0),  # offset added on to the Ainc values for all the defining sections in this surface
        "xles": np.array([0.0, 0.0]),  # leading edge cordinate vector(x component)
        "yles": np.array([-5.0, 0.0]),  # leading edge cordinate vector(y component)
        "zles": np.array([0.0, 0.0]),  # leading edge cordinate vector(z component)
        "chords": np.array([1.0, 1.0]),  # chord length vector
        "aincs": np.ones(2),#np.array([0.0, 0.0]),  # incidence angle vector

        # Geometry: Cross Sections
        # "xfminmax": np.array([[0.0, 1.0], [0.0, 1.0]]),  # airfoil x/c limits
        # NACA
        # 'naca' : np.array(['2412','2412']), # 4-digit NACA airfoil
        # Coordinates
        # 'xasec': np.array([[0., 1.], [0., 1.]]), # the x coordinate aifoil section
        # 'casec': np.array([[0., 0.], [0., 0.]]), # camber line at xasec
        # 'tasec': np.array([[0., 0.], [0., 0.]]), # thickness at xasec
        # 'xuasec': np.array([[0., 0.], [0., 0.]]), # airfoil upper surface x-coords (alternative to specifying camber line)
        # 'xlasec': np.array([[0., 0.], [0., 0.]]),  # airfoil lower surface x-coords (alternative to specifying camber line)
        # 'zuasec': np.array([[0., 0.], [0., 0.]]),  # airfoil upper surface z-coords (alternative to specifying camber line)
        # 'zlasec': np.array([[0., 0.], [0., 0.]]),  # airfoil lower surface z-coords (alternative to specifying camber line)
        # Airfoil Files
        'afiles': np.array(['airfoils/ag40d.dat','airfoils/ag40d.dat']), # airfoil file names

        # Paneling
        "nchordwise": np.int32(10),  # number of chordwise horseshoe vortice s placed on the surface
        "cspace": np.float64(0.0),  # chordwise vortex spacing parameter
        "nspan": np.int32(6),  # number of spanwise horseshoe vortices placed on the entire surface
        "sspace": np.float64(0.0),  # spanwise vortex spacing parameter for entire surface
        # "nspans": np.array([5, 5], dtype=np.int32),  # number of spanwise elements vector
        # "sspaces": np.array([3.0, 3.0], dtype=np.float64),  # spanwise spacing vector (for each section)
        "use surface spacing": np.int32(
            1
        ),  # surface spacing set under the surface heeading (known as LSURFSPACING in AVL)

        # Control Surfaces
        "control_assignments": {
            "flap" : {"assignment":np.array([0, 1],dtype=np.int32),
                      "xhinged": 0.8, # x/c location of hinge
                      "vhinged": np.zeros(3), # vector giving hinge axis about which surface rotates
                      "gaind": 1.0, # control surface gain
                      "refld": 1.0  # control surface reflection, sign of deflection for duplicated surface
                      }
        },

        # Design Variables (AVL)
        "design_var_assignments": {
            "des" : {"assignment":np.array([0, 1],dtype=np.int32),
                     "gaing":1.0}
        },

    }
}


geom_mesh = {
    "title": "Aircraft",
    "mach": np.float64(0.0),
    "iysym": np.int32(0),
    "izsym": np.int32(0),
    "zsym": np.float64(0.0),
    "Sref": np.float64(10.0),
    "Cref": np.float64(1.0),
    "Bref": np.float64(10.0),
    "XYZref": np.array([0.25, 0, 0],dtype=np.float64),
    "CDp": np.float64(0.0),
    "surfaces": surf,
    # Global Control and DV info
    "dname": ["flap"],  # Name of control input for each corresonding index
    "gname": ["des"],  # Name of design var for each corresonding index
}

geom_avl = {
    "title": "Aircraft",
    "mach": np.float64(0.0),
    "iysym": np.int32(0),
    "izsym": np.int32(0),
    "zsym": np.float64(0.0),
    "Sref": np.float64(10.0),
    "Cref": np.float64(1.0),
    "Bref": np.float64(10.0),
    "XYZref": np.array([0.25, 0, 0],dtype=np.float64),
    "CDp": np.float64(0.0),
    "surfaces": surf_avl,
    # Global Control and DV info
    "dname": ["flap"],  # Name of control input for each corresonding index
    "gname": ["des"],  # Name of design var for each corresonding index
}

keys_forces = ["CL", "CD"]

class TestMesh(unittest.TestCase):
    def setUp(self):
        self.ovl_mesh = OVLSolver(input_dict=geom_mesh)
        self.ovl_avl = OVLSolver(input_dict=geom_avl)

    def test_forces(self):
        self.ovl_mesh.set_variable("alpha", 2.0)
        self.ovl_avl.set_variable("alpha", 2.0)

        self.ovl_mesh.execute_run()
        self.ovl_avl.execute_run()

        forces_mesh = self.ovl_mesh.get_total_forces()
        forces_avl = self.ovl_avl.get_total_forces()

        for key in keys_forces:
            np.testing.assert_allclose(
                    forces_mesh[key],
                    forces_avl[key],
                    rtol=1e-8,
                )

    def test_control_surfaces(self):

        self.ovl_mesh.set_variable("alpha", 0.0)
        self.ovl_avl.set_variable("alpha", 0.0)

        self.ovl_mesh.set_control_deflection("flap", 2.0)
        self.ovl_avl.set_control_deflection("flap", 2.0)

        self.ovl_mesh.execute_run()
        self.ovl_avl.execute_run()

        forces_mesh = self.ovl_mesh.get_total_forces()
        forces_avl = self.ovl_avl.get_total_forces()

        for key in keys_forces:
            np.testing.assert_allclose(
                    forces_mesh[key],
                    forces_avl[key],
                    rtol=1e-8,
                )

if __name__ == "__main__":
    unittest.main()
