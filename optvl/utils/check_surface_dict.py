import warnings


def pre_check_input_dict(inputDict:dict):
    """
    This routine performs some verifications on a user's input diciontary to OptVL.
    It checks to see if any unsupported keys are in the inputs dictionary and the surface and body subdictionaries and issues a warning if any are detected.
    Additionally, we verify that the users has not specified redundant inputs for airfoil cross-section specification and y-symmetry sepcification and raise errors if any are detected.

    This routine does NOT type check inputs as that is handled in the load routine itself.

    NOTE: There are other redundant specifications in AVL where specifying one option will override what is specified in another.
    This overriding behavior is standard in AVL but in the future OptVL may check for these redundancies and raise warnings or errors.

    List non-checked redundancies:
    1. nspan and sspace override nspans and sspaces only if 'use surface spacing' is True
    2. 'cdcl' overrides 'cdclsec'


    Parameters
    ----------
    surface : dict
        User-defined OptVL input dict
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
            "Xref",
            "Yref",
            "Zref",
            "CDp",
            "surfaces",
            "bodies",
    ]

    # NOTE: make sure this is consistent to the documentation  page
    keys_implemented_surface = [
            #General
            "num_sections",
            'component', # logical surface component index (for grouping interacting surfaces, see AVL manual)
            'yduplicate', # surface is duplicated over the ysymm plane
            'nowake', # specifies that this surface is to NOT shed a wake, so that its strips will not have their Kutta conditions imposed
            'noalbe', # specifies that this surface is unaffected by freestream direction changes specified by the alpha,beta angles and p,q,r rotation rates
            'noload', # specifies that the force and moment on this surface is to NOT be included in the overall forces and moments of the configuration
            'clcdsec', # profile-drag CD(CL) function for each section in this surface
            'cdcl', # profile-drag CD(CL) function for all sections in this surface, overrides clcdsec.
            'claf', # CL alpha (dCL/da) scaling factor per section
            # Geometry
            'scale', # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
            'translate', # offset added on to all X,Y,Z values in this surface
            'angle', # offset added on to the Ainc values for all the defining sections in this surface
            'xles', # leading edge cordinate vector(x component)
            'yles', # leading edge cordinate vector(y component)
            'zles', # leading edge cordinate vector(z component)
            'chords', # chord length vector
            'aincs', # incidence angle vector
            # Geometry: Cross Sections
            'xfminmax', # airfoil x/c limits
            # NACA
            'naca', # 4-digit NACA airfoil
            # Manually specify airfoil coordinates in dictionary
            'airfoils',
            # Manual airfoil geometry
            'xasec', # the x coordinate aifoil section
            'casec', # camber line at xasec
            'tasec', # thickness at xasec
            'xuasec', # airfoil upper surface x-coords (alternative to specifying camber line)
            'xlasec',  # airfoil lower surface x-coords (alternative to specifying camber line)
            'zuasec',  # airfoil upper surface z-coords (alternative to specifying camber line)
            'zlasec',  # airfoil lower surface z-coords (alternative to specifying camber line)
            # Airfoil Files
            'afiles', # airfoil file names
            # Paneling
            'nchordwise', # number of chordwise horseshoe vortice s placed on the surface
            'cspace', # chordwise vortex spacing parameter
            'nspan', # number of spanwise horseshoe vortices placed on the entire surface
            'sspace', # spanwise vortex spacing parameter for entire surface
            'nspans', # number of spanwise elements vector, overriden by nspans
            'sspaces', # spanwise spacing vector (for each section), overriden by sspace
            'use surface spacing', # surface spacing set under the surface heeading (known as LSURFSPACING in AVL)
            # Control Surfaces
            # "dname" # IMPLEMENT THIS
            'icontd', # control variable index
            'xhinged', # x/c location of hinge
            'vhinged',  # vector giving hinge axis about which surface rotates
            'gaind', # control surface gain
            'refld', # control surface reflection, sign of deflection for duplicated surface
            'idestd', # design variable index
            'gaing', # desgin variable gain
    ]

    # NOTE: make sure this is consistent to the documentation  page
    keys_implemented_body = [
    ]

    multi_section_keys = [
        'nspans', # number of spanwise elements vector, overriden by nspans
        'sspaces', # spanwise spacing vector (for each section), overriden by sspace
        'clcdsec', # profile-drag CD(CL) function for each section in this surface
        'claf', # CL alpha (dCL/da) scaling factor per section
        # Geometry: Cross Sections
        # NACA
        'naca',
        # Coordinates
        'xasec', # the x coordinate aifoil section
        'casec', # camber line at xasec
        'tasec', # thickness at xasec
        'xuasec', # airfoil upper surface x-coords (alternative to specifying camber line)
        'xlasec',  # airfoil lower surface x-coords (alternative to specifying camber line)
        'zuasec',  # airfoil upper surface z-coords (alternative to specifying camber line)
        'zlasec',  # airfoil lower surface z-coords (alternative to specifying camber line)
        # Airfoil Files
        'afiles', # airfoil file names
        'xfminmax', # airfoil x/c limits
        # Paneling
        'nspans',
        'sspaces',
    ]

    # NOTE: make sure this is consistent to the documentation  page
    # Options used to specify airfoil sections for surfaces
    af_load_ops = ["naca", "airfoils", "afiles"]
    manual_af_override = ['xasec','casec','tasec','xuasec','xlasec','zuasec','zlasec']



    for key in inputDict.keys():
        # Check for keys not implemented
        if key not in keys_implemented_general:
            warnings.warn(
                "Key `{}` in input dict is (likely) not supported in OptVL and will be ignored".format(key),
                category=RuntimeWarning,
                stacklevel=2,
            )

    if ("surfaces" in inputDict.keys()):
        if (len(inputDict["surfaces"])>0):
            for surface in inputDict["surfaces"].keys():
                for key in inputDict["surfaces"][surface].keys():
                    # Check to verify if redundant y-symmetry specification are not made
                    if ("ydupl" in key) and ("iysym" in inputDict.keys()):
                        if (inputDict["surfaces"][surface]["yduplicate"] == 0.0) and (inputDict["iysym"] != 0):
                            raise RuntimeError(f"ERROR: Redundant y-symmetry specifications in surface {surface} \nIYSYM /= 0 \nYDUPLICATE  0.0. \nCan use one or the other, but not both!")

                    # Basically checks to see that at most only one of the options in af_load_ops or one of the options in manual_af_override is selected
                    if sum(bool(g) for g in ((af_load_ops & inputDict["surfaces"][surface].keys())|{any(k in manual_af_override for k in inputDict["surfaces"][surface].keys())})) > 1:
                        raise RuntimeError("More than one airfoil section specification detected in input dictionary!\n"
                        "Select only a single approach for specifying airfoil sections!")

                    # Check the surface input size is a 2D array with second dim equal to num_sections
                    if key in multi_section_keys:
                        if (key == "clcdsec") and (inputDict["surfaces"][surface][key].ndim != 3):
                            raise ValueError(f"Key {key} is of dimension {inputDict['surfaces'][surface][key].ndim}, expected 3!")
                        if (key != "clcdsec") and  inputDict["surfaces"][surface][key].ndim != 2:
                            raise ValueError(f"Key {key} is of dimension {inputDict['surfaces'][surface][key].ndim}, expected 2!")

                        if inputDict["surfaces"][surface][key].shape[1] != inputDict["surfaces"][surface]["num_sections"]:
                            raise ValueError(f"Key {key} does not have entries corresponding to each section!s")

                    # Check for keys not implemented
                    if key not in keys_implemented_surface:
                        warnings.warn(
                            "Key `{}` in surface dict {} is (likely) not supported in OptVL and will be ignored".format(key,surface),
                            category=RuntimeWarning,
                            stacklevel=2,
                        )

    if ("bodies" in inputDict.keys()):
        if (len(inputDict["bodies"])>0):
            for body in inputDict["bodies"].keys():
                for key in inputDict["bodies"][body].keys():
                    # Check to verify if redundant y-symmetry specification are not made
                    if ("ydupl" in key) and ("iysym" in inputDict.keys()):
                        if (inputDict["bodies"][body]["yduplicate"] == 0.0) and (inputDict["iysym"] != 0):
                            raise RuntimeError(f"ERROR: Redundant y-symmetry specifications in body {body} \nIYSYM /= 0 \nYDUPLICATE  0.0. \nCan use one or the other, but not both!")

                    # Basically checks to see that at most only one of the options in af_load_ops or one of the options in manual_af_override is selected
                    if sum(bool(g) for g in ((af_load_ops & inputDict["bodies"][body].keys())|{any(k in manual_af_override for k in inputDict["bodies"][body].keys())})) > 1:
                        raise RuntimeError("More than one airfoil section specification detected in input dictionary!\n"
                        "Select only a single approach for specifying airfoil sections!")

                    # Check for keys not implemented
                    if key not in keys_implemented_body:
                        warnings.warn(
                            "Key `{}` in body dict {} is (likely) not supported in OptVL and will be ignored".format(key,body),
                            category=RuntimeWarning,
                            stacklevel=2,
                        )

