"""
This module contains the routines, for pre-checking OptVL input dictionaries for errors.
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import warnings

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np

# =============================================================================
# Extension modules
# =============================================================================

def scalar_to_strip_vec(given_val,num_secs):
    """Converts a scalar into a numpy array of length matching
    the number of sections. If a numpy array of length not matching
    the number of sections is input then an exception is thrown.

    Args:
        given_val: Input value
        num_secs: Number of sections

    Returns:
        np.ndarray of length number of sections
    """
    # check if we input a scalar
    if not isinstance(given_val, np.ndarray):
        return given_val*np.ones(num_secs)
    elif isinstance(given_val,np.ndarray):
        if given_val.shape[0] != num_secs:
            raise ValueError("The length of a given surface/body input must either be a scalar or match the number of sections!")
        return given_val
    else:
        return given_val


def pre_check_input_dict(input_dict: dict):
    """
    This routine performs some verifications on a user's input diciontary to OptVL.
    It veries that geometry and control surfaces are specified in the correct format.
    Also checks to see if any unsupported keys are in the inputs dictionary and the surface and body subdictionaries and issues a warning if any are detected.
    
    This routine does NOT type check inputs as that is handled in the load routine itself.

    NOTE: There are other redundant specifications in AVL where specifying one option will override what is specified in another.
    This overriding behavior is standard in AVL but in the future OptVL may check for these redundancies and raise warnings or errors.

    List non-checked redundancies:
    1. nspan and sspace override nspans and sspaces only if 'use surface spacing' is True
    2. 'cdcl' overrides 'cdclsec'


    Args:
        input_dict : dict

    Returns:
        input_dict : dict
    """

    # NOTE: make sure this is consistent to the documentation  page
    keys_implemented_general = [
        "title",
        "mach",
        "iysym",
        "izsym",
        "zsym",
        "Sref",
        "Cref",
        "Bref",
        "XYZref",
        "CDp",
        "surfaces",
        "bodies",
        "dname",
        "gname",
    ]

    # NOTE: make sure this is consistent to the documentation  page
    keys_implemented_surface = [
        # General
        "num_sections",
        "num_controls",
        "num_design_vars",
        "component",  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate",  # surface is duplicated over the ysymm plane
        "wake",  # specifies that this surface is to NOT shed a wake, so that its strips will not have their Kutta conditions imposed
        "albe",  # specifies that this surface is unaffected by freestream direction changes specified by the alpha,beta angles and p,q,r rotation rates
        "load",  # specifies that the force and moment on this surface is to NOT be included in the overall forces and moments of the configuration
        "clcdsec",  # profile-drag CD(CL) function for each section in this surface
        "cdcl",  # profile-drag CD(CL) function for all sections in this surface, overrides clcdsec.
        "claf",  # CL alpha (dCL/da) scaling factor per section
        # Geometry
        "scale",  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate",  # offset added on to all X,Y,Z values in this surface
        "angle",  # offset added on to the Ainc values for all the defining sections in this surface
        "xles",  # leading edge cordinate vector(x component)
        "yles",  # leading edge cordinate vector(y component)
        "zles",  # leading edge cordinate vector(z component)
        "chords",  # chord length vector
        "aincs",  # incidence angle vector
        # Geometry: Cross Sections
        "xfminmax",  # airfoil x/c limits
        # NACA
        "naca",  # 4-digit NACA airfoil
        # Manually specify airfoil coordinates in dictionary
        "airfoils",
        # Manual airfoil geometry
        "xasec",  # the x coordinate aifoil section
        "casec",  # camber line at xasec
        "tasec",  # thickness at xasec
        "xuasec",  # airfoil upper surface x-coords (alternative to specifying camber line)
        "xlasec",  # airfoil lower surface x-coords (alternative to specifying camber line)
        "zuasec",  # airfoil upper surface z-coords (alternative to specifying camber line)
        "zlasec",  # airfoil lower surface z-coords (alternative to specifying camber line)
        # Airfoil Files
        "afiles",  # airfoil file names
        # Paneling
        "nchordwise",  # number of chordwise horseshoe vortice s placed on the surface
        "cspace",  # chordwise vortex spacing parameter
        "nspan",  # number of spanwise horseshoe vortices placed on the entire surface
        "sspace",  # spanwise vortex spacing parameter for entire surface
        "nspans",  # number of spanwise elements vector, overriden by nspans
        "sspaces",  # spanwise spacing vector (for each section), overriden by sspace
        "use surface spacing",  # surface spacing set under the surface heeading (known as LSURFSPACING in AVL)
        # Geometery: Mesh
        "mesh",
        "flatten_mesh",
        # Control Surfaces
        "control_assignments",
        "icontd",  # control variable index
        "xhinged",  # x/c location of hinge
        "vhinged",  # vector giving hinge axis about which surface rotates
        "gaind",  # control surface gain
        "refld",  # control surface reflection, sign of deflection for duplicated surface
        # Design Variables
        "design_var_assignments",
        "idestd",  # design variable index
        "gaing",  # desgin variable gain
    ]

    # NOTE: make sure this is consistent to the documentation  page
    keys_implemented_body = [
        "nvb",  # number of sources
        "bspace",  # source spacing
        "yduplicate",  # duplicate body over y-axis
        "scale",  # scaling factors applied to all x,y,z coordinates
        "translate",  # offset added on to all X,Y,Z values in this body
        "body_oml",
        "bfile",
    ]

    multi_section_keys = [
        "nspans",  # number of spanwise elements vector, overriden by nspans
        "sspaces",  # spanwise spacing vector (for each section), overriden by sspace
        "clcdsec",  # profile-drag CD(CL) function for each section in this surface
        "claf",  # CL alpha (dCL/da) scaling factor per section
        # Geometry
        "xles",  # leading edge cordinate vector(x component)
        "yles",  # leading edge cordinate vector(y component)
        "zles",  # leading edge cordinate vector(z component)
        "chords",  # chord length vector
        "aincs",  # incidence angle vector
        # Geometry: Cross Sections
        # NACA
        "naca",
        # Coordinates
        "xasec",  # the x coordinate aifoil section
        "casec",  # camber line at xasec
        "tasec",  # thickness at xasec
        "xuasec",  # airfoil upper surface x-coords (alternative to specifying camber line)
        "xlasec",  # airfoil lower surface x-coords (alternative to specifying camber line)
        "zuasec",  # airfoil upper surface z-coords (alternative to specifying camber line)
        "zlasec",  # airfoil lower surface z-coords (alternative to specifying camber line)
        # Airfoil Files
        "afiles",  # airfoil file names
        "xfminmax",  # airfoil x/c limits
        # Paneling
        "nspans",
        "sspaces",
    ]

    control_keys = [
        "icontd",  # control variable index
        "xhinged",  # x/c location of hinge
        "vhinged",  # vector giving hinge axis about which surface rotates
        "gaind",  # control surface gain
        "refld",  # control surface reflection, sign of deflection for duplicated surface
    ]

    design_var_keys =[
        "idestd",
        "gaing",  # design variable surface gain
    ]


    dim_2_keys = [
        "clcdsec",
        "xfminmax",
        "xasec",
        "casec",
        "tasec",
        "xuasec",
        "xlasec",
        "zuasec",
        "zlasec",
    ]

    # NOTE: make sure this is consistent to the documentation  page
    # Options used to specify airfoil sections for surfaces
    airfoil_spec_keys = ["naca", "airfoils", "afiles", "xasec"]

    for key in input_dict.keys():

        # Check if the user provided negative reference values
        if key in ["Bref", "Sref", "Cref"]:
            if input_dict[key] < 0.0:
                raise ValueError(f"Reference value {key} cannot be negative!")

        # Correct incorrect symmetry plane defs with warning
        if key in ["iysym", "izsym"]:
            if input_dict[key] not in [-1,0,1]:
                warnings.warn(
                    f"OptVL WARNING - Option {key} needs to be -1, 0, or 1!\n "
                    f"Correcting by setting based on sign to {np.sign(input_dict[key])}.\n",
                    stacklevel=2,
                )
                input_dict[key] = np.sign(input_dict[key])

        # Check for keys not implemented
        if key not in keys_implemented_general:
            warnings.warn(
                "Key `{}` in input dict is (likely) not supported in OptVL and will be ignored".format(key),
                category=RuntimeWarning,
                stacklevel=2,
            )

    if "surfaces" in input_dict.keys():
        if len(input_dict["surfaces"]) > 0:
            for surface in input_dict["surfaces"].keys():

                # Check if we are directly providing a mesh and set the strips as "sections" so that the maps setup correctly
                if "mesh" in input_dict["surfaces"][surface].keys():
                    # First check if the mesh is a valid numpy array shape
                    if len(input_dict["surfaces"][surface]["mesh"].shape) != 3:
                        raise ValueError("The provided mesh must be a numpy array of size (nx,ny,3)")
                    # If we are using a mesh then set number of sections equal to number of strip for the purposes of intialization
                    input_dict["surfaces"][surface]["num_sections"] = input_dict["surfaces"][surface]["mesh"].shape[1]

                # Verify at least two section
                if input_dict["surfaces"][surface]["num_sections"] < 2:
                    raise RuntimeError("Must have at least two sections per surface!")

                # Read and process the controls dictionary
                if "control_assignments" in input_dict["surfaces"][surface]:
                    num_controls_per_sec = np.zeros(input_dict["surfaces"][surface]["num_sections"],dtype=np.int32)

                    for control in input_dict["surfaces"][surface]["control_assignments"]:
                        if control not in input_dict["dname"]:
                            raise ValueError(f"Control {control}, in surface {surface} not defined in dname!")

                        # built the control data lists if needed 
                        if "icontd" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["icontd"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]
                        if "xhinged" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["xhinged"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]
                        if "vhinged" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["vhinged"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]
                        if "gaind" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["gaind"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]
                        if "refld" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["refld"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]

                        # Add one to the number of controls defined for each section
                        sec_assign = input_dict["surfaces"][surface]["control_assignments"][control]["assignment"]
                        num_controls_per_sec[sec_assign] += 1
                        # assign data to sections
                        for idx_sec in input_dict["surfaces"][surface]["control_assignments"][control]["assignment"]:
                            input_dict["surfaces"][surface]["icontd"][idx_sec].append(input_dict["dname"].index(control)) # Add control index to icontd for each section
                            input_dict["surfaces"][surface]["xhinged"][idx_sec].append(input_dict["surfaces"][surface]["control_assignments"][control]["xhinged"]) # Add hinge line position
                            input_dict["surfaces"][surface]["vhinged"][idx_sec].append(input_dict["surfaces"][surface]["control_assignments"][control]["vhinged"]) # Add hinge vector position
                            input_dict["surfaces"][surface]["gaind"][idx_sec].append(input_dict["surfaces"][surface]["control_assignments"][control]["gaind"]) # Add gain information
                            input_dict["surfaces"][surface]["refld"][idx_sec].append(input_dict["surfaces"][surface]["control_assignments"][control]["refld"]) # Add reflection information

                    # set the control numbers per section
                    input_dict["surfaces"][surface]["num_controls"] = num_controls_per_sec
                elif "num_controls" not in input_dict["surfaces"][surface]:
                    # Otherwise if we are not manually specifying controls then zero out the num_controls array
                    input_dict["surfaces"][surface]["num_controls"] = np.zeros(input_dict["surfaces"][surface]["num_sections"],dtype=np.int32)

                 # Read and process the design variables dictionary
                if "design_var_assignments" in input_dict["surfaces"][surface]:
                    num_design_vars_per_sec = np.zeros(input_dict["surfaces"][surface]["num_sections"],dtype=np.int32)

                    for design_var in input_dict["surfaces"][surface]["design_var_assignments"]:
                        if design_var not in input_dict["gname"]:
                            raise ValueError(f"Design Variable {design_var}, in surface {surface} not defined in gname!")

                        # built the control data lists if needed 
                        if "idestd" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["idestd"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]
                        if "gaing" not in input_dict["surfaces"][surface]:
                            input_dict["surfaces"][surface]["gaing"] = [[] for _ in range(input_dict["surfaces"][surface]["num_sections"])]


                        # Add one to the number of controls defined for each section
                        sec_assign = input_dict["surfaces"][surface]["design_var_assignments"][design_var]["assignment"]
                        num_design_vars_per_sec[sec_assign] += 1
                        # assign data to sections
                        for idx_sec in input_dict["surfaces"][surface]["design_var_assignments"][design_var]["assignment"]:
                            input_dict["surfaces"][surface]["idestd"][idx_sec].append(input_dict["gname"].index(design_var)) # Add design var index to idestd for each section
                            input_dict["surfaces"][surface]["gaing"][idx_sec].append(input_dict["surfaces"][surface]["design_var_assignments"][design_var]["gaing"]) # Add gain information

                    # set the control numbers per section
                    input_dict["surfaces"][surface]["num_design_vars"] = num_design_vars_per_sec
                elif "num_design_vars" not in input_dict["surfaces"][surface]:
                    # Otherwise if we are not manually specifying controls then zero out the num_design_vars array
                    input_dict["surfaces"][surface]["num_design_vars"] = np.zeros(input_dict["surfaces"][surface]["num_sections"],dtype=np.int32)


                #Checks to see that at most only one of the options in af_load_ops or one of the options in manual_af_override is selected
                if len(airfoil_spec_keys & input_dict["surfaces"][surface].keys()) > 1:
                    raise RuntimeError(
                        "More than one airfoil section specification detected in input dictionary!\n"
                        "Select only a single approach for specifying airfoil sections!")

                # Process all keys
                for key in input_dict["surfaces"][surface].keys():

                    # Check to verify if redundant y-symmetry specification are not made
                    if ("ydupl" in key) and ("iysym" in input_dict.keys()):
                        if (input_dict["surfaces"][surface]["yduplicate"] == 0.0) and (input_dict["iysym"] != 0):
                            raise RuntimeError(
                                f"ERROR: Redundant y-symmetry specifications in surface {surface} \nIYSYM /= 0 \nYDUPLICATE  0.0. \nCan use one or the other, but not both!"
                            )

                    # Verify that keys that need items specified for every strip/section have value specified for all strip/section or have a scalar/single vector that can be duplicated
                    if key in multi_section_keys:

                        if (key in dim_2_keys):
                            if not (isinstance(input_dict["surfaces"][surface][key],np.ndarray)):
                                raise ValueError(f"Input for {key} must be a single dim 1 numpy array or a dim 2 numpy array with each vector along axis 0 corresponding to strip/section")

                            # If the user provides a single dim 1 vector stack it num_sections times
                            if (input_dict["surfaces"][surface][key].ndim == 1):
                                input_dict["surfaces"][surface][key] = np.tile(input_dict["surfaces"][surface][key],(input_dict["surfaces"][surface]["num_sections"],1))
                            # Otherwise make sure we have entries for each seciton
                            elif (input_dict["surfaces"][surface][key].ndim == 2):
                                if (input_dict["surfaces"][surface][key].shape[0] != input_dict["surfaces"][surface]["num_sections"]):
                                    raise ValueError(
                                        f"Key {key} only has {input_dict['surfaces'][surface][key].shape[0]}, expected {input_dict['surfaces'][surface]['num_sections']}!"
                                    )
                            else:
                                raise ValueError(
                                f"Key {key} is of dimension {input_dict['surfaces'][surface][key].ndim}, expected 1 or 2!"
                            )
                        else:
    
                            # If the user provides a scalar or string expand it out for all sections
                            if isinstance(input_dict["surfaces"][surface][key],(int,float,np.int32,np.float64,str)):
                                    input_dict["surfaces"][surface][key] = np.tile(input_dict["surfaces"][surface][key],(input_dict["surfaces"][surface]["num_sections"]))
                            elif input_dict["surfaces"][surface][key].ndim > 1:
                                raise ValueError(
                                    f"Key {key} is of dimension {input_dict['surfaces'][surface][key].ndim}, expected 1!"
                                )

                    # Check for keys not implemented
                    if key not in keys_implemented_surface:
                        warnings.warn(
                            "Key `{}` in surface dict {} is (likely) not supported in OptVL and will be ignored".format(
                                key, surface
                            ),
                            category=RuntimeWarning,
                            stacklevel=2,
                        )

                    # Check if controls defined correctly
                    if key in control_keys:
                        for j in range(input_dict["surfaces"][surface]["num_sections"]):
                            for _ in range(input_dict["surfaces"][surface]["num_controls"][j]):
                                if (
                                    len(input_dict["surfaces"][surface][key][j])
                                    != input_dict["surfaces"][surface]["num_controls"][j]
                                ):
                                    raise ValueError(
                                        f"Key {key} does not have entries corresponding to each control for this section!"
                                    )


                    # Check if dvs defined correctly
                    if key in design_var_keys:
                        for j in range(input_dict["surfaces"][surface]["num_sections"]):
                            for _ in range(input_dict["surfaces"][surface]["num_design_vars"][j]):
                                if (
                                    len(input_dict["surfaces"][surface][key][j])
                                    != input_dict["surfaces"][surface]["num_design_vars"][j]
                                ):
                                    raise ValueError(
                                        f"Key {key} does not have entries corresponding to each design var for this section!"
                                    )


    else:
        # Add dummy entry if surfaces are not defined
        input_dict["surfaces"] = {}

    if "bodies" in input_dict.keys():
        if len(input_dict["bodies"]) > 0:
            for body in input_dict["bodies"].keys():
                # Check that only one body oml input is selected
                if ("body_oml" in input_dict["bodies"][body].keys()) and ("bfile" in input_dict["bodies"][body].keys()):
                    raise RuntimeError("Select only one body oml definition!")
                elif ("body_oml" not in input_dict["bodies"][body].keys()) and (
                    "bfile" not in input_dict["bodies"][body].keys()
                ):
                    raise RuntimeError("Must define a oml for a body!")

                for key in input_dict["bodies"][body].keys():
                    # Check to verify if redundant y-symmetry specification are not made
                    if ("ydupl" in key) and ("iysym" in input_dict.keys()):
                        if (input_dict["bodies"][body]["yduplicate"] == 0.0) and (input_dict["iysym"] != 0):
                            raise RuntimeError(
                                f"ERROR: Redundant y-symmetry specifications in body {body} \nIYSYM /= 0 \nYDUPLICATE  0.0. \nCan use one or the other, but not both!"
                            )

                    # Check if user tried to use body sections
                    if key == "num_sections":
                        raise RuntimeError(
                            "Body sections are a cut feature from AVL and are hence not support in OptVL."
                        )

                    # Check for keys not implemented
                    if key not in keys_implemented_body:
                        warnings.warn(
                            "Key `{}` in body dict {} is (likely) not supported in OptVL and will be ignored".format(
                                key, body
                            ),
                            category=RuntimeWarning,
                            stacklevel=2,
                        )
    else:
        # Add dummy entry if bodies are not defined
        input_dict["bodies"] = {}

    return input_dict
