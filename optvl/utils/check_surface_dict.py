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
        # Control Surfaces
        # "dname" # IMPLEMENT THIS
        "icontd",  # control variable index
        "xhinged",  # x/c location of hinge
        "vhinged",  # vector giving hinge axis about which surface rotates
        "gaind",  # control surface gain
        "refld",  # control surface reflection, sign of deflection for duplicated surface
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
    total_global_control = 0
    total_global_design_var = 0
    if "surfaces" in input_dict.keys():
        if len(input_dict["surfaces"]) > 0:
            for surface in input_dict["surfaces"].keys():

                # Verify at least two section
                if input_dict["surfaces"][surface]["num_sections"] < 2:
                    raise RuntimeError("Must have at least two sections per surface!")

                # if no controls are specified then fill it in with 0s
                if "num_controls" not in input_dict["surfaces"][surface].keys():
                    input_dict["surfaces"][surface]["num_controls"] = np.zeros(input_dict["surfaces"][surface]["num_sections"],dtype=np.int32)

                 # if no dvs are specified then fill it in with 0s
                if "num_design_vars" not in input_dict["surfaces"][surface].keys():
                    input_dict["surfaces"][surface]["num_design_vars"] = np.zeros(input_dict["surfaces"][surface]["num_sections"],dtype=np.int32)
                
                #Checks to see that at most only one of the options in af_load_ops or one of the options in manual_af_override is selected
                if len(airfoil_spec_keys & input_dict["surfaces"][surface].keys()) > 1:
                    raise RuntimeError(
                        "More than one airfoil section specification detected in input dictionary!\n"
                        "Select only a single approach for specifying airfoil sections!")

                
                for key in input_dict["surfaces"][surface].keys():

                    # Check to verify if redundant y-symmetry specification are not made
                    if ("ydupl" in key) and ("iysym" in input_dict.keys()):
                        if (input_dict["surfaces"][surface]["yduplicate"] == 0.0) and (input_dict["iysym"] != 0):
                            raise RuntimeError(
                                f"ERROR: Redundant y-symmetry specifications in surface {surface} \nIYSYM /= 0 \nYDUPLICATE  0.0. \nCan use one or the other, but not both!"
                            )

                    # Check the surface input size is a 2D array with second dim equal to num_sections
                    if key in multi_section_keys:
                        if (key in dim_2_keys) and (input_dict["surfaces"][surface][key].ndim != 2):
                            raise ValueError(
                                f"Key {key} is of dimension {input_dict['surfaces'][surface][key].ndim}, expected 2!"
                            )
                        if (key not in dim_2_keys) and input_dict["surfaces"][surface][key].ndim != 1:
                            raise ValueError(
                                f"Key {key} is of dimension {input_dict['surfaces'][surface][key].ndim}, expected 1!"
                            )

                        if (
                            input_dict["surfaces"][surface][key].shape[0]
                            != input_dict["surfaces"][surface]["num_sections"]
                        ):
                            raise ValueError(f"Key {key} does not have entries corresponding to each section!s")

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
                                    input_dict["surfaces"][surface][key][j].shape[0]
                                    != input_dict["surfaces"][surface]["num_controls"][j]
                                ):
                                    raise ValueError(
                                        f"Key {key} does not have entries corresponding to each control for this section!"
                                    )

                    # Accumulate icont max
                    if "icontd" in input_dict["surfaces"][surface].keys():
                        arr = input_dict["surfaces"][surface]["icontd"]
                        vals = [a.max() + 1 for a in arr if a.size > 0]
                        total_global_control = max(vals) if vals else None
                        # total_global_control = np.max(input_dict["surfaces"][surface]["icontd"])+1

                    # Check if dvs defined correctly
                    if key in control_keys:
                        for j in range(input_dict["surfaces"][surface]["num_sections"]):
                            for _ in range(input_dict["surfaces"][surface]["num_design_vars"][j]):
                                if (
                                    input_dict["surfaces"][surface][key][j].shape[0]
                                    != input_dict["surfaces"][surface]["num_design_vars"][j]
                                ):
                                    raise ValueError(
                                        f"Key {key} does not have entries corresponding to each design var for this section!"
                                    )

                    # Accumulate idestd max
                    if "idestd" in input_dict["surfaces"][surface].keys():
                        arr = input_dict["surfaces"][surface]["idestd"]
                        vals = [a.max() + 1 for a in arr if a.size > 0]
                        total_global_design_var = max(vals) if vals else None
                        # total_global_design_var = np.max(input_dict["surfaces"][surface]["idestd"])+1

            if "icontd" in input_dict["surfaces"][surface].keys():
                if len(input_dict["dname"]) != (total_global_control):
                    raise ValueError(
                        "Number of unique control names does not match the number of unique controls defined!"
                    )

            if "idestd" in input_dict["surfaces"][surface].keys():
                if len(input_dict["gname"]) != (total_global_design_var):
                    raise ValueError(
                        "Number of unique design vars does not match the number of unique controls defined!"
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
