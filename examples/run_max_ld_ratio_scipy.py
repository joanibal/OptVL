"""
Script to perform L/D maximization for an aircraft using OptVL and Scipy.
It optimizes Angle of Attack (alpha) and Aspect Ratio (AR).
"""

import numpy as np
from scipy.optimize import minimize
# Assuming OptVL is imported as 'optvl' and the solver class is 'OVLSolver'
# from optvl import OVLSolver  # Or the correct import path and class name
# For now, defining a placeholder if the actual class is not available in this environment
class OVLSolverPlaceholder:
    def __init__(self, geo_file, debug=False):
        self.geo_file = geo_file
        self.debug = debug
        self.parameters = {}
        self.current_totals = {}
        self.current_sensitivities = {}
        print(f"Placeholder OVLSolver initialized with {geo_file}, debug={debug}")

    def load_geometry(self, geo_file=None):
        if geo_file:
            self.geo_file = geo_file
        print(f"Placeholder: Loading geometry {self.geo_file}...")
        # Simulate loading some default values for Sref, b, Cref for placeholder
        self.current_totals = {'Sref': 20.0, 'b': 10.0, 'Cref': 2.0, 'CL': 0.1, 'CD': 0.02}

    def set_parameter(self, *args):
        # Simplified: store param path and value
        # e.g., set_parameter("Mach", 0.1) -> self.parameters["Mach"] = 0.1
        # e.g., set_parameter("SURFACE", 1, "SECTION", 1, "SPAN", 5.0)
        #       -> self.parameters["SURFACE_1_SECTION_1_SPAN"] = 5.0
        key = args[0]
        if len(args) > 2:
            key = "_".join(map(str, args[:-1]))
        self.parameters[key] = args[-1]
        # print(f"Placeholder: Setting parameter {key} = {args[-1]}")

    def execute_run(self):
        print("Placeholder: Executing run...")
        # Simulate some behavior based on alpha and AR if they are set
        alpha = self.parameters.get("alpha", 0.0)
        # AR is not directly used by AVL, but span/chord are.
        # This is a very simplified placeholder effect.
        cl_factor = (alpha / 5.0) + 0.1 # crude alpha effect
        cd_factor = 0.02 + (alpha / 5.0)**2 * 0.01 # crude alpha effect
        
        # Simulate effect of SPAN/CHORD changes (which are driven by AR)
        # This is highly simplified and not aerodynamically accurate.
        span_param = self.parameters.get("SURFACE_1_SECTION_1_SPAN", B_INITIAL / 2.0 if 'B_INITIAL' in globals() else 5.0)
        chord_param = self.parameters.get("SURFACE_1_SECTION_1_CHORD", C_INITIAL if 'C_INITIAL' in globals() else 2.0)
        
        # Effective AR from params for simulation
        sim_ar = (2 * span_param)**2 / ( (2*span_param) * chord_param) if chord_param > 1e-6 else 1.0
        
        self.current_totals['CL'] = 0.1 + cl_factor + (sim_ar / 8.0 - 1.0) * 0.2
        self.current_totals['CD'] = 0.02 + cd_factor + (8.0 / sim_ar - 1.0) * 0.01
        # Ensure Sref, b, Cref are present from initial load
        self.current_totals.setdefault('Sref', SREF_INITIAL if 'SREF_INITIAL' in globals() else 20.0)
        self.current_totals.setdefault('b', B_INITIAL if 'B_INITIAL' in globals() else 10.0)
        self.current_totals.setdefault('Cref', C_INITIAL if 'C_INITIAL' in globals() else 2.0)
        return self.current_totals


    def get_total_forces(self):
        print("Placeholder: Getting total forces...")
        return self.current_totals

    def execute_run_sensitivities(self, var_list):
        print(f"Placeholder: Executing run sensitivities for {var_list}...")
        # Simulate some sensitivities
        # These would depend on the current state (alpha, AR derived params)
        alpha_val = self.parameters.get("alpha", 0.0)
        span_param_val = self.parameters.get("SURFACE_1_SECTION_1_SPAN", B_INITIAL / 2.0 if 'B_INITIAL' in globals() else 5.0)
        chord_param_val = self.parameters.get("SURFACE_1_SECTION_1_CHORD", C_INITIAL if 'C_INITIAL' in globals() else 2.0)

        self.current_sensitivities = {
            'CL': {
                'alpha': 0.1, # dCL/dalpha per degree (typical 2*pi/180 * rad_to_deg_factor ~ 0.1)
                'SURFACE_1_SECTION_1_SPAN': 0.05, # dCL/d(span_param)
                'SURFACE_1_SECTION_1_CHORD': 0.15  # dCL/d(chord_param)
            },
            'CD': {
                'alpha': 0.005, # dCD/dalpha
                'SURFACE_1_SECTION_1_SPAN': -0.002, # dCD/d(span_param)
                'SURFACE_1_SECTION_1_CHORD': 0.01   # dCD/d(chord_param)
            }
        }
        return self.current_sensitivities

# Global storage for initial geometric properties
SREF_INITIAL = None
B_INITIAL = None
C_INITIAL = None # Mean Geometric Chord

# Initialize OptVL Solver
# Replace OVLSolverPlaceholder with actual OptVL solver class
# For example: from optvl.solver import OVLSolver
try:
    from optvl.solver import OVLSolver # Attempt to import the actual solver
    print("Successfully imported OVLSolver from optvl.solver")
except ImportError:
    print("Failed to import OVLSolver from optvl.solver, using placeholder.")
    OVLSolver = OVLSolverPlaceholder # Use placeholder if import fails

ovl_solver = OVLSolver(geo_file="aircraft.avl", debug=False)

def setup_initial_conditions():
    """
    Loads geometry, sets initial Mach, runs a baseline analysis,
    and stores initial geometric properties.
    """
    global SREF_INITIAL, B_INITIAL, C_INITIAL
    ovl_solver.load_geometry("aircraft.avl") # Ensure correct path if not in examples/
    ovl_solver.set_parameter("Mach", 0.1)
    ovl_solver.set_parameter("alpha", 0.0) # Nominal alpha for baseline run
    
    # Execute run to populate initial aerodynamic and geometric data
    # In a real scenario, OVLSolver.load_geometry might return these,
    # or they are fetched after load_geometry.
    # If execute_run() is needed to get Sref, b, Cref, it implies they might be calculated values.
    initial_run_data = ovl_solver.execute_run() # Assuming this updates get_total_forces() source
    initial_totals = ovl_solver.get_total_forces()

    if not initial_totals:
        raise ValueError("Failed to get initial total forces from OVLSolver.")

    SREF_INITIAL = initial_totals.get('Sref')
    B_INITIAL = initial_totals.get('b') # Initial total span
    C_INITIAL = initial_totals.get('Cref') # Initial reference chord (often mean geometric chord)

    if SREF_INITIAL is None or B_INITIAL is None or C_INITIAL is None:
        raise ValueError(
            f"Could not retrieve initial geometric properties: "
            f"Sref={SREF_INITIAL}, b={B_INITIAL}, Cref={C_INITIAL}. "
            "Check OVLSolver's get_total_forces() output format or load_geometry behavior."
        )
    print(f"Initial conditions set: Sref={SREF_INITIAL}, Span={B_INITIAL}, Chord_ref={C_INITIAL}")

def apply_design_variables(x):
    """Applies design variables [alpha_deg, aspect_ratio] to the ovl_solver."""
    alpha_deg, current_AR = x
    
    ovl_solver.set_parameter("alpha", alpha_deg)

    # Aspect Ratio Modification (maintaining SREF_INITIAL)
    # new_span is the total span for the new AR
    new_span = np.sqrt(current_AR * SREF_INITIAL)
    # new_chord is the new mean geometric chord for the new AR and SREF_INITIAL
    new_chord = SREF_INITIAL / new_span 

    # Apply to OptVL:
    # These parameter names and the strategy depend heavily on:
    # 1. The structure of "aircraft.avl" (e.g., is it a simple rectangular wing, tapered, etc.?)
    # 2. How OVLSolver's `set_parameter` maps to AVL's geometry modification capabilities.
    # Assumption: Modifying SURFACE 1, SECTION 1's SPAN and CHORD parameters.
    # - 'SPAN' for a section in AVL typically refers to the semi-span of that segment if YLEAD is not used.
    #   If SECTION 1 defines the root and there's a SECTION 2 for the tip,
    #   then SECTION 2's YLE (Y leading edge) would be new_span / 2.0.
    #   And SECTION 1 YLE would be 0.
    # - 'CHORD' for a section is its local chord.
    # For a simple rectangular wing where SURFACE 1 SECTION 1 defines the whole half-wing:
    #   SECTION 1 'SPAN' parameter would be new_span / 2.0
    #   SECTION 1 'CHORD' parameter would be new_chord
    # This is a simplification. A robust implementation needs detailed knowledge of the AVL file
    # and how OptVL manipulates it (e.g., scaling all sections, modifying specific section YLEs/Chords).
    
    # --- BEGIN COMMENT: Geometry Modification Assumptions ---
    # The following `set_parameter` calls for SPAN and CHORD are based on these assumptions:
    # 1. `aircraft.avl` has a primary lifting surface designated as 'SURFACE 1'.
    # 2. This surface has at least one 'SECTION 1'.
    # 3. The 'SPAN' parameter for 'SECTION 1' controls its semi-span (half of the total span B_INITIAL for a simple wing).
    #    Thus, we set it to `new_span / 2.0`.
    # 4. The 'CHORD' parameter for 'SECTION 1' controls its chord. For a rectangular wing, this would be `new_chord`.
    #    If the wing is tapered, this would likely be the root chord, and other sections would need
    #    to be scaled accordingly, or OptVL handles this scaling internally based on a root chord change.
    # This part is highly dependent on the specific `aircraft.avl` file's structure and the OptVL API.
    # A more general approach might involve scaling all X, Y, Z coordinates of wing sections.
    # --- END COMMENT: Geometry Modification Assumptions ---
    ovl_solver.set_parameter('SURFACE', 1, 'SECTION', 1, 'SPAN', new_span / 2.0)
    ovl_solver.set_parameter('SURFACE', 1, 'SECTION', 1, 'CHORD', new_chord)


def objective_function(x):
    """
    Calculates -L/D ratio for a given set of design variables x = [alpha_deg, aspect_ratio].
    """
    apply_design_variables(x)
    
    ovl_solver.execute_run()
    totals = ovl_solver.get_total_forces()

    CL = totals.get('CL', 0.0)
    CD = totals.get('CD', 1e6) # Default to a large number if CD is not found

    if CD < 1e-8: # Avoid division by zero or very small CD
        L_D_ratio = -1e6 # Penalize
    else:
        L_D_ratio = CL / CD
    
    print(f"Objective: x={x}, L/D={L_D_ratio:.4f} (CL={CL:.4f}, CD={CD:.5f})")
    return -L_D_ratio # Minimize -L/D to maximize L/D

def objective_gradient(x):
    """
    Calculates the gradient of -L/D ratio w.r.t. x = [alpha_deg, aspect_ratio].
    """
    alpha_deg, current_AR = x

    # Ensure solver state is correct for sensitivity analysis
    apply_design_variables(x)
    ovl_solver.execute_run() # Ensure AVL is in the correct state FOR sensitivities

    # Get sensitivities of CL and CD w.r.t. available OptVL variables
    # The sensitivity keys must match what OptVL provides, which in turn should relate
    # to the parameters that can be set by `set_parameter`.
    sens = ovl_solver.execute_run_sensitivities(['CL', 'CD'])

    if not sens or 'CL' not in sens or 'CD' not in sens:
        print(f"Warning: Could not retrieve sensitivities for x={x}. Returning zero gradient.")
        return np.array([0.0, 0.0])

    # 1. Gradients w.r.t. alpha
    dCL_dalpha = sens['CL'].get('alpha', 0.0)
    dCD_dalpha = sens['CD'].get('alpha', 0.0)

    # 2. Gradients w.r.t. Aspect Ratio (AR) - Implemented via Chain Rule
    #    OptVL provides sensitivities w.r.t. direct geometry parameters like 'SPAN' or 'CHORD'
    #    of sections. We need to map these to AR.

    #    Let `span_param` be the value set by `ovl_solver.set_parameter('SURFACE', 1, 'SECTION', 1, 'SPAN', ...)`,
    #    which is `new_span / 2.0`.
    #    Let `chord_param` be the value set by `ovl_solver.set_parameter('SURFACE', 1, 'SECTION', 1, 'CHORD', ...)`,
    #    which is `new_chord`.

    #    Define the conceptual names for these parameters as keys in the sensitivity dictionary:
    #    (These must match exactly how OVLSolver returns sensitivity keys)
    span_param_key = 'SURFACE_1_SECTION_1_SPAN' # e.g., sensitivity w.r.t. value passed to set_parameter
    chord_param_key = 'SURFACE_1_SECTION_1_CHORD'# e.g., sensitivity w.r.t. value passed to set_parameter

    dCL_dSPAN_param = sens['CL'].get(span_param_key, 0.0)
    dCD_dSPAN_param = sens['CD'].get(span_param_key, 0.0)
    dCL_dCHORD_param = sens['CL'].get(chord_param_key, 0.0)
    dCD_dCHORD_param = sens['CD'].get(chord_param_key, 0.0)

    # Chain Rule Derivatives:
    # current_AR = x[1]
    # new_span (b) = sqrt(current_AR * SREF_INITIAL)
    # new_chord (c) = SREF_INITIAL / new_span = sqrt(SREF_INITIAL / current_AR)
    # span_param (let's call it s_p) = new_span / 2.0
    # chord_param (let's call it c_p) = new_chord
    
    # db_dAR = 0.5 * (current_AR * SREF_INITIAL)**(-0.5) * SREF_INITIAL
    #        = 0.5 * SREF_INITIAL / new_span = 0.5 * new_chord
    db_dAR = 0.5 * new_chord # where new_chord is sqrt(SREF_INITIAL / current_AR) at current_AR

    # dc_dAR = d(sqrt(SREF_INITIAL / current_AR))/dAR
    #        = sqrt(SREF_INITIAL) * (-0.5) * current_AR**(-1.5)
    #        = -0.5 * sqrt(SREF_INITIAL / current_AR**3)
    #        = -0.5 * (sqrt(SREF_INITIAL / current_AR) / current_AR) = -0.5 * new_chord / current_AR
    dc_dAR = -0.5 * new_chord / current_AR
    
    # d(span_param)/dAR = d(new_span/2.0)/dAR = 0.5 * db_dAR
    dSPAN_param_dAR = 0.5 * db_dAR # = 0.25 * new_chord

    # d(chord_param)/dAR = d(new_chord)/dAR = dc_dAR
    dCHORD_param_dAR = dc_dAR # = -0.5 * new_chord / current_AR

    # --- BEGIN COMMENT: Chain Rule for d/dAR ---
    # The calculation of dCL_dAR and dCD_dAR uses the chain rule:
    # dF/dAR = (dF/dSPAN_param * dSPAN_param/dAR) + (dF/dCHORD_param * dCHORD_param/dAR)
    # where F is CL or CD.
    # - dSPAN_param/dAR and dCHORD_param/dAR are the derivatives of the actual values
    #   passed to `set_parameter` (i.e., new_span/2 and new_chord) with respect to AR.
    # - This assumes that `OVLSolver` provides sensitivities `dCL/dSPAN_param` and `dCL/dCHORD_param`
    #   (and similarly for CD) corresponding to the `span_param_key` and `chord_param_key`.
    # - If OptVL offers direct sensitivities dCL/dAR and dCD/dAR (e.g., if AR itself can be
    #   a "variable" input to a geometry engine within OptVL), that would be far simpler and preferred.
    #   The current implementation assumes sensitivities are w.r.t. lower-level geometry inputs.
    # --- END COMMENT: Chain Rule for d/dAR ---
    dCL_dAR = dCL_dSPAN_param * dSPAN_param_dAR + dCL_dCHORD_param * dCHORD_param_dAR
    dCD_dAR = dCD_dSPAN_param * dSPAN_param_dAR + dCD_dCHORD_param * dCHORD_param_dAR
    
    # Gradient of the objective function F = -CL/CD
    # dF/dvar = - (CD * dCL/dvar - CL * dCD/dvar) / CD^2
    #         =   (CL * dCD/dvar - CD * dCL/dvar) / CD^2
    totals = ovl_solver.get_total_forces() # Get current CL, CD
    CL = totals.get('CL', 0.0)
    CD = totals.get('CD', 1e-8) # Use small positive if CD is zero or not found, to avoid div by zero

    if abs(CD) < 1e-8:
        obj_grad_alpha = 0.0
        obj_grad_AR = 0.0
    else:
        obj_grad_alpha = (CL * dCD_dalpha - CD * dCL_dalpha) / (CD**2)
        obj_grad_AR    = (CL * dCD_dAR    - CD * dCL_dAR)    / (CD**2)
    
    gradient = np.array([obj_grad_alpha, obj_grad_AR])
    print(f"Gradient: x={x}, grad={gradient} (dCL/da={dCL_dalpha:.3f}, dCD/da={dCD_dalpha:.3f}, dCL/dAR={dCL_dAR:.3f}, dCD/dAR={dCD_dAR:.3f})")
    return gradient

# --- Main Execution Block ---
if __name__ == "__main__":
    print("Starting L/D Maximization Script...")
    try:
        setup_initial_conditions()
    except ValueError as e:
        print(f"Error during setup: {e}")
        print("Exiting script.")
        exit()
    except Exception as e: # Catch any other OVLSolver or file errors
        print(f"An unexpected error occurred during setup: {e}")
        print("Exiting script.")
        exit()


    # Initial guess: [alpha_deg, aspect_ratio]
    # Use B_INITIAL and SREF_INITIAL to calculate initial AR
    initial_AR = (B_INITIAL**2) / SREF_INITIAL if SREF_INITIAL > 1e-6 else 8.0
    x0 = np.array([2.0, initial_AR]) 
    print(f"Initial Guess (x0): Alpha={x0[0]} deg, AR={x0[1]:.2f}")

    # Bounds for design variables: (alpha_bounds, AR_bounds)
    bounds = ((-5.0, 10.0), (4.0, 12.0)) 
    print(f"Bounds: Alpha={bounds[0]}, AR={bounds[1]}")

    print("\nStarting optimization with Scipy SLSQP...")
    result = minimize(
        objective_function, 
        x0, 
        method='SLSQP', 
        jac=objective_gradient, 
        bounds=bounds, 
        options={'disp': True, 'ftol': 1e-7, 'maxiter': 50}, # ftol for objective change tolerance
        tol=1e-7 # tol for overall termination criteria (depends on solver)
    )

    print("\n--- Optimization Results ---")
    if result.success:
        optimal_vars = result.x
        max_L_D_ratio = -result.fun # Objective was -L/D
        print(f"Optimization successful!")
        print(f"  Optimal Alpha:        {optimal_vars[0]:.4f} degrees")
        print(f"  Optimal Aspect Ratio: {optimal_vars[1]:.4f}")
        print(f"  Maximum L/D Ratio:    {max_L_D_ratio:.4f}")

        # Verification run with optimal variables
        print("\n--- Verification Run with Optimal Variables ---")
        apply_design_variables(optimal_vars)
        ovl_solver.execute_run()
        final_totals = ovl_solver.get_total_forces()
        final_CL = final_totals.get('CL', 0.0)
        final_CD = final_totals.get('CD', 1e6)
        final_L_D = final_CL / final_CD if final_CD > 1e-8 else -1e6
        
        print(f"  Verified CL: {final_CL:.4f}")
        print(f"  Verified CD: {final_CD:.5f}")
        print(f"  Verified L/D Ratio: {final_L_D:.4f}")
    else:
        print("Optimization failed.")
        print(f"  Message: {result.message}")
        print(f"  Status: {result.status}")

    print("\nScript finished.")

```
