"""
Lift-Constrained Drag Minimization of a Wing using OptVL and Scipy.

This script demonstrates how to minimize the drag coefficient (CD) of a wing
subject to a constraint on the lift coefficient (CL).

Tutorial Steps:
1. Define the aerodynamic performance evaluation function using OptVL.
2. Create wrapper functions for Scipy's optimizer.
3. Set up and run the optimization.
4. Print and verify the results.

Assumptions about `examples/rectangle.avl`:
- It contains a primary lifting surface (assumed to be SURFACE 1).
- This surface has at least two sections (assumed SECTION 1 for root, SECTION 2 for tip for twist).
- OptVL can modify parameters like 'alpha', surface 'SPAN', 'CHORD', and section 'TWIST'.
- OptVL can return 'SREF', 'SPAN', 'CHORD' for a given surface after loading a file.
"""

import numpy as np
from scipy.optimize import minimize
import optvl # Assuming optvl is the import name for the OptVL library

# --- Constants ---
AVL_FILE_PATH = "examples/rectangle.avl" # Path to the AVL geometry file
TARGET_CL = 0.5
MACH_NUMBER = 0.1 # Using a low Mach number as AVL is primarily for incompressible flow
# Reynolds number can also be set if needed, e.g., optvl_instance.set_reynolds_number(1e6)

# Store initial geometric properties globally or pass them around carefully
initial_geom_properties = {}

def get_initial_geometry_properties(optvl_instance, avl_file):
    """
    Loads the AVL file and extracts initial geometric properties of the main wing.
    Stores them in the global `initial_geom_properties` dictionary.

    Args:
        optvl_instance: An initialized OptVL object.
        avl_file (str): Path to the .avl file.

    Returns:
        bool: True if properties were successfully extracted, False otherwise.
    """
    try:
        optvl_instance.load_geometry(avl_file)
        # Assumptions:
        # - The main wing is SURFACE 1.
        # - OptVL provides methods like `get_value` to fetch these parameters.
        #   The exact names ('SREF', 'SPAN', 'CHORD') might vary based on OptVL's API.
        sref = optvl_instance.get_value("SURFACE", 1, "SREF")
        span = optvl_instance.get_value("SURFACE", 1, "SPAN")
        chord = optvl_instance.get_value("SURFACE", 1, "CHORD") # Assuming this is mean geometric chord

        if sref is None or span is None or chord is None:
            print("Error: Could not retrieve initial Sref, Span, or Chord from OptVL.")
            print(f"  Sref: {sref}, Span: {span}, Chord: {chord}")
            return False

        initial_geom_properties['sref'] = sref
        initial_geom_properties['span'] = span
        initial_geom_properties['chord'] = chord
        initial_geom_properties['ar'] = span**2 / sref # Calculate initial AR
        print(f"Initial properties: Sref={sref:.4f}, Span={span:.4f}, Chord={chord:.4f}, AR={initial_geom_properties['ar']:.4f}")
        return True
    except Exception as e:
        print(f"Error in get_initial_geometry_properties: {e}")
        return False

def modify_wing_geometry(optvl_instance, aspect_ratio, wing_twist_deg):
    """
    Modifies the wing geometry (aspect ratio and twist) in the OptVL instance.
    Aspect ratio is changed while keeping Sref constant.

    Args:
        optvl_instance: An OptVL object with geometry loaded.
        aspect_ratio (float): The target aspect ratio.
        wing_twist_deg (float): The linear twist at the wing tip (degrees). Root twist is 0.

    Returns:
        bool: True if modification was successful, False otherwise.
    """
    try:
        sref_initial = initial_geom_properties['sref']

        # Calculate new span and chord to achieve the target AR while preserving Sref
        new_span = np.sqrt(aspect_ratio * sref_initial)
        # new_chord = sref_initial / new_span # This would be the new mean geometric chord

        # Apply changes. The exact OptVL API calls might differ.
        # Assumption: OptVL allows direct modification of surface span.
        # If OptVL requires scaling factors, they would be:
        # span_scale = new_span / initial_geom_properties['span']
        # chord_scale = new_chord / initial_geom_properties['chord']
        # And then apply these scales to X, Y, Z coordinates of sections.
        # For simplicity, we assume direct span setting is possible.
        # If OptVL changes chord automatically when span and Sref are set, that's ideal.
        # Or, if we need to set chord as well:
        # optvl_instance.set_value("SURFACE", 1, "CHORD", new_chord)
        
        # This is a key part: How OptVL handles geometry modification.
        # Option 1: Set span, and OptVL adjusts chord to maintain Sref (if Sref is fixed by OptVL)
        # Option 2: Set span, and also set chord.
        # Option 3: Scale all sections. For a simple rectangle, scaling X and Y of sections.
        # Let's assume OptVL provides a way to set span and chord directly for a surface,
        # or that setting span and having Sref fixed (if that's an OptVL concept) is enough.
        # For this tutorial, we'll assume setting SPAN and CHORD is the OptVL way.
        # However, for a simple rectangular wing, changing span and keeping Sref constant means
        # all sectional chords must change proportionally.
        # `optvl_instance.set_value("SURFACE", 1, "CHORD", new_chord_value_for_all_sections)`
        # A more robust way for generic AVL files would be to scale section chords.
        
        # Let's assume we can set the overall surface span.
        # And we'll adjust the CHORD parameter for *all* sections of SURFACE 1.
        # This assumes rectangle.avl defines its chord via SECTION CHORD values.
        optvl_instance.set_value("SURFACE", 1, "SPAN", new_span)
        
        num_sections = optvl_instance.get_value("SURFACE", 1, "NUMSECTIONS") # Needs API support
        if num_sections is None: # Fallback if NUMSECTIONS is not available
            print("Warning: Could not get NUMSECTIONS. Assuming 2 sections for rectangle.avl.")
            num_sections = 2 # Typical for a simple rectangle.avl file (root and tip)

        new_section_chord = sref_initial / new_span # This is the new chord for all sections of a rect wing
        for i in range(1, num_sections + 1):
             optvl_instance.set_value("SURFACE", 1, "SECTION", i, "CHORD", new_section_chord)


        # Apply wing twist
        # Assuming linear twist: 0 at root (SECTION 1), `wing_twist_deg` at tip (SECTION N)
        # If there are intermediate sections, their twist would be interpolated by AVL.
        optvl_instance.set_value("SURFACE", 1, "SECTION", 1, "TWIST", 0.0)
        # Assuming the last section defines the tip for twist purposes.
        optvl_instance.set_value("SURFACE", 1, "SECTION", num_sections, "TWIST", wing_twist_deg)
        
        return True
    except Exception as e:
        print(f"Error in modify_wing_geometry: {e}")
        return False


def get_aero_performance(design_vars, optvl_instance_template):
    """
    Calculates aerodynamic performance (CD, CL) for a given set of design variables.

    Args:
        design_vars (list or np.array): [alpha_deg, aspect_ratio, wing_twist_deg].
        optvl_instance_template: An initialized OptVL object (used for cloning or re-init).

    Returns:
        tuple: (CD, CL). Returns (np.inf, -np.inf) or similar on error.
    """
    alpha_deg, aspect_ratio, wing_twist_deg = design_vars

    # Create a fresh OptVL instance or reload to ensure clean state
    # This depends on OptVL's behavior. For safety, reloading is good.
    current_optvl = optvl.OptVL() # Or clone optvl_instance_template if possible
    try:
        current_optvl.load_geometry(AVL_FILE_PATH)
        current_optvl.set_mach_number(MACH_NUMBER)
        # current_optvl.set_reynolds_number(REYNOLDS_NUMBER) # If desired

        # Set angle of attack
        current_optvl.set_parameter("alpha", alpha_deg)

        # Modify geometry
        if not modify_wing_geometry(current_optvl, aspect_ratio, wing_twist_deg):
            print("Error: Geometry modification failed.")
            return np.inf, -np.inf # Penalize if geometry modification fails

        # Execute AVL run
        results = current_optvl.execute_run()

        if results and 'totals' in results and 'CDtot' in results['totals'] and 'CLtot' in results['totals']:
            cd = results['totals']['CDtot']
            cl = results['totals']['CLtot']
            # print(f"Vars: a={alpha_deg:.2f}, AR={aspect_ratio:.2f}, Tw={wing_twist_deg:.2f} -> CD={cd:.5f}, CL={cl:.5f}")
            return cd, cl
        else:
            print("Error: OptVL execution failed or results format unexpected.")
            print(f"  Results: {results}")
            return np.inf, -np.inf # Penalize if results are not as expected

    except Exception as e:
        print(f"Exception in get_aero_performance for vars {design_vars}: {e}")
        return np.inf, -np.inf # Penalize heavily on any error
    finally:
        # Clean up OptVL instance if necessary (e.g., current_optvl.close())
        pass


def objective_function(design_vars, optvl_instance_template):
    """Scipy objective function: returns CD."""
    cd, _ = get_aero_performance(design_vars, optvl_instance_template)
    return cd

def constraint_function(design_vars, optvl_instance_template):
    """Scipy constraint function: returns CL - TARGET_CL (for equality constraint CL == TARGET_CL)."""
    _, cl = get_aero_performance(design_vars, optvl_instance_template)
    return cl - TARGET_CL

# --- Main Optimization Script ---
if __name__ == "__main__":
    print("Starting Lift-Constrained Drag Minimization Tutorial Script...")

    # Initialize OptVL
    # This instance is mainly for getting initial properties and as a template.
    base_optvl = optvl.OptVL()
    try:
        if not get_initial_geometry_properties(base_optvl, AVL_FILE_PATH):
            raise RuntimeError("Failed to get initial geometry properties. Exiting.")
    except Exception as e:
        print(f"Failed to initialize OptVL or load initial geometry: {e}")
        exit()

    # Design Variables: [alpha_deg, aspect_ratio, wing_twist_deg]
    # Initial guess (x0)
    x0 = [2.0, initial_geom_properties.get('ar', 8.0), 0.0] # Alpha (deg), AR, Wing Tip Twist (deg)

    # Bounds for design variables
    # Alpha: -5 to 10 deg
    # AR: 4 to 12 (example)
    # Twist: -5 to 2 deg (example)
    bounds = [(-5.0, 10.0), (4.0, 12.0), (-5.0, 2.0)]

    # Constraints for Scipy
    # We want CL = TARGET_CL, so constraint_function should return 0.
    constraints = [{'type': 'eq', 'fun': constraint_function, 'args': (base_optvl,)}]

    print(f"\nInitial Guess (x0): Alpha={x0[0]} deg, AR={x0[1]}, Twist={x0[2]} deg")
    print(f"Bounds: {bounds}")
    print(f"Target CL: {TARGET_CL}")
    print(f"Optimizer: SLSQP")

    # Run optimization
    print("\nStarting optimization...")
    try:
        result = minimize(objective_function, x0, args=(base_optvl,),
                          method='SLSQP',
                          bounds=bounds,
                          constraints=constraints,
                          options={'disp': True, 'ftol': 1e-6, 'maxiter': 100})

        print("\n--- Optimization Results ---")
        if result.success:
            print("Optimization successful!")
            optimal_vars = result.x
            optimal_cd = result.fun
            # Verify CL with optimal variables
            # For verification, it's best to get CL from a direct call, not from the constraint function's last run
            _, optimal_cl_verify = get_aero_performance(optimal_vars, base_optvl)

            print(f"  Optimal Alpha:          {optimal_vars[0]:.4f} degrees")
            print(f"  Optimal Aspect Ratio:   {optimal_vars[1]:.4f}")
            print(f"  Optimal Wing Twist:     {optimal_vars[2]:.4f} degrees")
            print(f"  Minimized CD:           {optimal_cd:.6f}")
            print(f"  Achieved CL (verify):   {optimal_cl_verify:.6f} (Target: {TARGET_CL})")

            # Further verification with a completely fresh instance
            print("\n--- Verification with Fresh OptVL Instance ---")
            verify_optvl = optvl.OptVL()
            final_cd, final_cl = get_aero_performance(optimal_vars, verify_optvl)
            print(f"  CD from fresh run: {final_cd:.6f}")
            print(f"  CL from fresh run: {final_cl:.6f}")
            if abs(final_cl - TARGET_CL) > 1e-3: # Looser tolerance for verification print
                 print(f"Warning: Final CL ({final_cl:.6f}) deviates significantly from target ({TARGET_CL}).")

        else:
            print("Optimization failed.")
            print(f"  Message: {result.message}")
            print(f"  Status: {result.status}")

    except Exception as e:
        print(f"An error occurred during optimization: {e}")

    finally:
        # Clean up base_optvl if needed (e.g., base_optvl.close())
        print("\nScript finished.")

"""
Example `examples/rectangle.avl` structure assumed:

SimpleRectangularWing
#Mach
0.1
#IYsym   IZsym   IXsym
0       0       0
#Sref   Cref    Bref
20.0    2.0     10.0  <-- These are reference values, OptVL might use them or calculate from geometry
#Xref   Yref    Zref
0.0     0.0     0.0
#CDp
0.0100
#
# WING (SURFACE 1)
SURFACE
MainWing
#Nchordwise  Cspace   Nspanwise   Sspace
12           1.0      20          1.0
#
# ROOT SECTION (SECTION 1 of SURFACE 1)
SECTION
#Xle    Yle    Zle    Chord   Ainc  Nspanwise Sspace
0.0     0.0     0.0    2.0     0.0   10        1.0  <-- Chord for this section
AFILEN
naca0012.dat
#
# TIP SECTION (SECTION 2 of SURFACE 1)
SECTION
#Xle    Yle    Zle    Chord   Ainc  Nspanwise Sspace
0.0     5.0     0.0    2.0     0.0   0         1.0  <-- Chord for this section, Yle defines half-span
AFILEN
naca0012.dat

Key aspects for geometry modification:
- Sref, Span, Chord are defined for the whole aircraft. OptVL must allow reading these.
- Individual sections have Xle, Yle, Zle, Chord, Ainc (twist).
- To change Aspect Ratio while keeping Sref:
    - new_span = sqrt(AR_new * Sref_initial)
    - new_section_chord = Sref_initial / new_span (for a rectangular wing)
    - Modify Yle of tip section to new_span / 2.
    - Modify Chord of all sections to new_section_chord.
- To change Twist:
    - Modify Ainc of sections. For linear twist from 0 at root to T at tip:
        - Root section (Yle=0): Ainc = 0
        - Tip section (Yle=new_span/2): Ainc = wing_twist_deg
- The script assumes `optvl_instance.set_value("SURFACE", 1, "SPAN", new_span)` correctly adjusts
  the underlying geometry (like section Yle values) and
  `optvl_instance.set_value("SURFACE", 1, "SECTION", i, "CHORD", new_section_chord)`
  updates section chords. If OptVL needs explicit Xle, Yle, Zle, Chord settings for each
  section, the `modify_wing_geometry` function would be more complex.
  The current implementation of `modify_wing_geometry` using `SPAN` and `CHORD` setters
  is a higher-level abstraction that depends on specific OptVL capabilities.
  A more robust implementation for generic AVL files would involve getting number of sections,
  then iterating through sections to update their YLE (for span) and CHORD.
  `optvl_instance.set_value("SURFACE", 1, "SECTION", <tip_sec_idx>, "YLE", new_span / 2.0)`
"""
