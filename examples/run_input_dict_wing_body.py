import numpy as np
from optvl import OVLSolver

# Example: Specify a wing body configuration using the input dictionary features and write the equivalent AVL input file

# Specify the surfaces: in this case only a wing
surf = {
    "Wing": {
        # General
        "num_sections": np.int32(2), # number of sections in surface
        "num_controls": np.array([1, 0], dtype=np.int32), # number of control surfaces assocaited with each section (in order)
        "num_design_vars": np.array([1, 0], dtype=np.int32), # number of AVL design variables assocaited with each section (in order)
        "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        "wake": np.int32(
            1
        ),  # specifies that this surface is to shed a wake, so that its strips will not have their Kutta conditions imposed (default true)
        "albe": np.int32(
            1
        ),  # specifies that this surface is affected by freestream direction changes specified by the alpha,beta angles and p,q,r rotation rates (default true)
        "load": np.int32(
            1
        ),  # specifies that the force and moment on this surface is to be included in the overall forces and moments of the configuration (default true)
        "clcdsec": np.array(
            [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        ),  # profile-drag CD(CL) function for each section in this surface (default 0.0)
        "cdcl": np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),  # profile-drag CD(CL) function for all sections in this surface, overrides clcdsec. (default 0.0)
        "claf": np.array([1.0, 1.0]),  # CL alpha (dCL/da) scaling factor for each section in order (default 1.0)
        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords are also scaled by Xscale)
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
        # CHOOSE ONLY ONE APPROACH
        "xfminmax": np.array([[0.0, 1.0], [0.0, 1.0]]),  # airfoil x/c limits
        # NACA
        # 'naca' : np.array(['2412','2412']), # 4-digit NACA airfoil
        # Airfoil
        #'airfoils' : np.array([[0., 1.],[0., 0.]]), # Specify airfoil coordinates for each section
        # Raw Geometry
        # 'xasec': np.array([[0., 1.], [0., 1.]]), # the x coordinate aifoil section
        # 'casec': np.array([[0., 0.], [0., 0.]]), # camber line at xasec
        # 'tasec': np.array([[0., 0.], [0., 0.]]), # thickness at xasec
        # 'xuasec': np.array([[0., 0.], [0., 0.]]), # airfoil upper surface x-coords
        # 'xlasec': np.array([[0., 0.], [0., 0.]]),  # airfoil lower surface x-coords
        # 'zuasec': np.array([[0., 0.], [0., 0.]]),  # airfoil upper surface z-coords
        # 'zlasec': np.array([[0., 0.], [0., 0.]]),  # airfoil lower surface z-coords
        # Airfoil Files
        'afiles': np.array(['airfoils/ag40d.dat','airfoils/ag40d.dat']), # airfoil file names
        # Paneling
        "nchordwise": np.int32(5),  # number of chordwise horseshoe vortice s placed on the surface
        "cspace": np.float64(0.0),  # chordwise vortex spacing parameter
        "nspan": np.int32(2),  # number of spanwise horseshoe vortices placed on the entire surface
        "sspace": np.float64(0.0),  # spanwise vortex spacing parameter for entire surface
        "nspans": np.array([5, 5], dtype=np.int32),  # number of spanwise elements vector
        "sspaces": np.array([3.0, 3.0], dtype=np.float64),  # spanwise spacing vector (for each section)
        "use surface spacing": np.int32(
            1
        ),  # surface spacing set for the entire surface (known as LSURFSPACING in AVL)
        # Control Surfaces
        "icontd": [np.array([0],dtype=np.int32), np.array([],dtype=np.int32)],  # control variable indicies assocaited with each section in order
        "xhinged": [np.array([0.8]), np.array([],dtype=np.float64)],  # x/c location of hinge for each control index associated with each section
        "vhinged": [np.array([[0.0, 0.0, 0.0]]), np.array([],dtype=np.float64)],  # vector giving hinge axis about which surface rotates for each control index associated with each section
        "gaind": [np.array([1.0]), np.array([],dtype=np.float64)],  # control surface gain for each control index associated with each section
        "refld": [np.array([1.0]), np.array([],dtype=np.float64)],  # control surface reflection, sign of deflection for duplicated surface for each control index associated with each section
        # Design Variables (AVL)
        "idestd": [np.array([0],dtype=np.int32), np.array([],dtype=np.int32)],  # design variable indicies assocaited with each section in order
        "gaing": [np.array([1.0]), np.array([],dtype=np.float64)],  # desgin variable gain for each control index associated with each section
    }
}


body = {
    "Fuselage": {
        # General
        # 'yduplicate': np.float64(0), # body is duplicated over the ysymm plane
        # Geometry
        "scale": np.array(
            [[1.0, 1.0, 1.0]]
        ),  # scaling factors applied to all x,y,z coordinates (chords areal so scaled by Xscale)
        "translate": np.array([[0.0, 0.0, 0.0]]),  # offset added on to all X,Y,Z values in this surface
        # Geometry: OML
        # CHOOSE ONLY ONE APPROACH
        # OML Files
        "bfile": "fuseSupra.dat",  # body oml file name
        #"body_oml": np.zeros(50), # Specify body coordiantes directly
        # Discretization
        "nvb": np.int32(10),  # number of source-line nodes
        "bspace": np.float64(0.0),  # lengthwise node spacing parameter
    }
}


inputDict = {
    "title": "Aircraft", # Aicraft name (MUST BE SET)
    "mach": np.float64(0.0), # Reference Mach number
    "iysym": np.int32(0), # y-symmetry settings
    "izsym": np.int32(0), # z-symmetry settings
    "zsymm": np.float64(0.0), # z-symmetry plane
    "Sref": np.float64(10.0), # Reference planform area
    "Cref": np.float64(1.0), # Reference chord area
    "Bref": np.float64(10.0), # Reference span length
    "XYZref": np.array([0.25, 0, 0],dtype=np.float64), # Reference x,y,z position
    "CDp": np.float64(0.0), # Reference profile drag adjustment
    "surfaces": surf, # dictionary of surface dictionaries
    "bodies": body, # dictionary of body dictionaries
    # Global Control and DV info
    "dname": np.array(["flap"]),  # Name of control input for each corresonding index
    "gname": np.array(["des"]),  # Name of design var for each corresonding index
}

ovl = OVLSolver(input_dict=inputDict, debug=True, timing=True)

# look at the geometry to see that everything is right
ovl.plot_geom()

# set the angle of attack
ovl.set_constraint("alpha", 2.0)

# set the deflection of the trailing edge flap
ovl.set_constraint("flap", 2.0)

# set mach number
ovl.set_parameter("Mach", 0.3)

# This is the method that acutally runs the analysis
ovl.execute_run()

# print data about the run
force_data = ovl.get_total_forces()
print(
    f'CL:{force_data["CL"]:10.6f}   CD:{force_data["CD"]:10.6f}   CM:{force_data["CM"]:10.6f}'
)

# lets look at the cp countours
ovl.plot_cp()

# Get the body parameters
body_data = ovl.get_body_params()

# Change some body parameters
body_data["Fuselage"]["nvb"] = np.int32(15)

# Set the changed parameters
ovl.set_body_params(body_data)

# Get the surface parameters
surf_data = ovl.get_surface_params(
    include_section_geom=True, include_paneling=True, include_con_surf=True, include_airfoils=True
)

# change the surface parameters
surf_data["Wing"]["nspan"] = np.int32(4)

# set the surface parameters
ovl.set_surface_params(surf_data)


# Run the analysis again
ovl.execute_run()

# lets look at the cp countours again
ovl.plot_cp()


# Get the whole input dict that would give us what we have in there
input_data = ovl.get_input_dict()

# Let's look at it
print(input_data)

# Lastly, let's write the classical AVL input file
ovl.write_geom_file("aircraft_wing_body_afile.avl")
