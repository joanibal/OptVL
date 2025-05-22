"""
Script to perform parameter sweeps (Angle of Attack and Aspect Ratio)
using OptVL and plot results with Matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt

# Global storage for initial geometric properties
SREF_INITIAL = None
B_INITIAL = None # Initial total span
C_INITIAL = None # Initial mean geometric chord (or reference chord)

# --- OVLSolver Placeholder (if real one is not available) ---
class OVLSolverPlaceholder:
    def __init__(self, geo_file, debug=False):
        self.geo_file = geo_file
        self.debug = debug
        self.parameters = {}
        self.current_totals = {'Sref': 20.0, 'b': 10.0, 'Cref': 2.0} # Default initial geo
        print(f"Placeholder OVLSolver initialized with {geo_file}, debug={debug}")

    def load_geometry(self, geo_file=None):
        if geo_file:
            self.geo_file = geo_file
        print(f"Placeholder: Loading geometry {self.geo_file}...")
        # Simulate loading some default values for Sref, b, Cref for placeholder
        # These will be updated by get_initial_geometry_properties
        self.current_totals = {'Sref': 20.0, 'b': 10.0, 'Cref': 2.0, 'CL': 0.1, 'CD': 0.02}


    def set_parameter(self, *args):
        key = args[0]
        if len(args) > 2: # e.g. ('SURFACE', 1, 'SECTION', 1, 'SPAN', val)
            key = "_".join(map(str, args[:-1]))
        self.parameters[key] = args[-1]
        # print(f"Placeholder: Setting parameter {key} = {args[-1]}")


    def execute_run(self):
        # print("Placeholder: Executing run...")
        alpha = self.parameters.get("alpha", 0.0)
        
        # Placeholder physics - crude simulation for plotting
        # Use initial AR if span/chord params not set yet, else calculate effective AR
        eff_AR = B_INITIAL**2 / SREF_INITIAL if SREF_INITIAL and B_INITIAL else 8.0
        
        # Try to get span and chord if they were set by AR modification logic
        # These keys must match those used in apply_design_variables
        span_param = self.parameters.get("SURFACE_1_SECTION_1_SPAN", B_INITIAL / 2.0 if B_INITIAL else 5.0)
        chord_param = self.parameters.get("SURFACE_1_SECTION_1_CHORD", C_INITIAL if C_INITIAL else 2.0)

        if SREF_INITIAL: # Ensure SREF_INITIAL is available
            calc_span = span_param * 2.0
            # For a rectangular wing, Sref = calc_span * chord_param
            # So, if chord_param is section chord, then effective AR:
            eff_AR = calc_span**2 / (calc_span * chord_param) if chord_param > 1e-6 else 1.0
        else: # Fallback if SREF_INITIAL not yet set (should not happen in normal flow)
             eff_AR = B_INITIAL**2 / SREF_INITIAL if SREF_INITIAL and B_INITIAL else 8.0

        # CL simulation: linear with alpha, slightly affected by AR
        cl_alpha_slope = 0.1 + (eff_AR - 8.0) * 0.005 # AR effect on lift slope
        self.current_totals['CL'] = cl_alpha_slope * alpha + (eff_AR / 8.0 - 1.0) * 0.05 # Base CL shift with AR

        # CD simulation: parabolic with alpha (induced drag like), base drag affected by AR (e.g. wetted area)
        cd0 = 0.02 - (eff_AR - 8.0) * 0.001 # AR effect on parasite drag
        self.current_totals['CD'] = cd0 + self.current_totals['CL']**2 / (np.pi * eff_AR * 0.8) # Induced drag like term
        
        # Ensure Sref, b, Cref are present from initial load or properties function
        self.current_totals.setdefault('Sref', SREF_INITIAL if SREF_INITIAL else 20.0)
        self.current_totals.setdefault('b', B_INITIAL if B_INITIAL else 10.0)
        self.current_totals.setdefault('Cref', C_INITIAL if C_INITIAL else 2.0)
        return self.current_totals

    def get_total_forces(self):
        # print("Placeholder: Getting total forces...")
        return self.current_totals

# Attempt to import actual OVLSolver
try:
    from optvl.solver import OVLSolver # Or actual import path
    print("Successfully imported OVLSolver.")
except ImportError:
    print("Failed to import OVLSolver. Using OVLSolverPlaceholder.")
    OVLSolver = OVLSolverPlaceholder

# --- Global OVLSolver instance ---
ovl_solver = OVLSolver(geo_file="aircraft.avl", debug=False)

def get_initial_geometry_properties():
    """
    Loads geometry, sets initial Mach, runs a baseline analysis,
    and stores initial geometric properties.
    """
    global SREF_INITIAL, B_INITIAL, C_INITIAL
    ovl_solver.load_geometry("aircraft.avl") # Path relative to where script is run
    ovl_solver.set_parameter("Mach", 0.1) # Default Mach for baseline
    ovl_solver.set_parameter("alpha", 0.0) # Nominal alpha for baseline run
    
    initial_run_data = ovl_solver.execute_run()
    initial_totals = ovl_solver.get_total_forces()

    if not initial_totals:
        raise ValueError("Failed to get initial total forces from OVLSolver for geometry properties.")

    SREF_INITIAL = initial_totals.get('Sref')
    B_INITIAL = initial_totals.get('b')       # Initial total span
    C_INITIAL = initial_totals.get('Cref')    # Initial reference chord (often mean geometric chord)

    if SREF_INITIAL is None or B_INITIAL is None or C_INITIAL is None:
        raise ValueError(
            f"Could not retrieve initial geometric properties: "
            f"Sref={SREF_INITIAL}, b={B_INITIAL}, Cref={C_INITIAL}. "
            "Check OVLSolver's get_total_forces() output or load_geometry behavior."
        )
    print(f"Initial properties set: Sref={SREF_INITIAL:.2f}, Span={B_INITIAL:.2f}, Chord_ref={C_INITIAL:.2f}")


def apply_aspect_ratio_modification(current_AR):
    """
    Modifies the wing geometry to achieve the target Aspect Ratio (AR)
    while maintaining the initial SREF_INITIAL.
    """
    if SREF_INITIAL is None:
        raise RuntimeError("Initial geometric properties not set. Call get_initial_geometry_properties() first.")

    new_span = np.sqrt(current_AR * SREF_INITIAL)
    new_chord = SREF_INITIAL / new_span # This is the new mean geometric chord

    # --- BEGIN COMMENT: Geometry Modification Assumptions ---
    # The following `set_parameter` calls for SPAN and CHORD are based on assumptions
    # about `aircraft.avl` structure and the OptVL API, similar to `run_max_ld_ratio_scipy.py`.
    # 1. `aircraft.avl` has a primary lifting surface as 'SURFACE 1'.
    # 2. This surface has at least one 'SECTION 1'.
    # 3. The 'SPAN' parameter for 'SECTION 1' controls its semi-span. We set it to `new_span / 2.0`.
    # 4. The 'CHORD' parameter for 'SECTION 1' controls its chord. For a rectangular wing, this is `new_chord`.
    #    If the wing is tapered, this would likely be the root chord, and OptVL/AVL would handle scaling.
    # This part is highly dependent on the specific `aircraft.avl` file and OptVL API.
    # --- END COMMENT: Geometry Modification Assumptions ---
    ovl_solver.set_parameter('SURFACE', 1, 'SECTION', 1, 'SPAN', new_span / 2.0)
    ovl_solver.set_parameter('SURFACE', 1, 'SECTION', 1, 'CHORD', new_chord)


# --- Main Script ---
if __name__ == "__main__":
    print("Starting Parameter Sweep Script...")

    try:
        get_initial_geometry_properties()
    except Exception as e:
        print(f"Error getting initial geometry properties: {e}")
        exit()

    # Define parameter ranges
    alpha_range = np.linspace(-2, 10, 13)  # 13 points from -2 to 10 degrees
    ar_range = np.array([6, 8, 10, 12])    # Discrete AR values

    print(f"\nSweeping Alpha from {alpha_range[0]} to {alpha_range[-1]} deg ({len(alpha_range)} points)")
    print(f"Sweeping Aspect Ratio for values: {ar_range}")

    results_data = []

    # Main loops for sweeping parameters
    for ar_val in ar_range:
        print(f"\nSetting Aspect Ratio = {ar_val}")
        try:
            apply_aspect_ratio_modification(ar_val)
        except Exception as e:
            print(f"  Error applying AR modification for AR={ar_val}: {e}. Skipping this AR.")
            continue

        for alpha_val in alpha_range:
            ovl_solver.set_parameter("alpha", alpha_val)
            
            try:
                ovl_solver.execute_run()
                totals = ovl_solver.get_total_forces()

                CL = totals.get('CL', np.nan)
                CD = totals.get('CD', np.nan)
                
                L_D_ratio = np.nan
                if CD > 1e-8 and not (np.isnan(CL) or np.isnan(CD)):
                    L_D_ratio = CL / CD
                
                results_data.append({
                    'alpha': alpha_val, 
                    'AR': ar_val, 
                    'CL': CL, 
                    'CD': CD, 
                    'L_D': L_D_ratio
                })
                print(f"  Alpha={alpha_val:5.2f} deg, AR={ar_val:4.1f} -> CL={CL:7.4f}, CD={CD:8.5f}, L/D={L_D_ratio:7.3f}")

            except Exception as e:
                print(f"  Error during OVLSolver run for Alpha={alpha_val}, AR={ar_val}: {e}")
                results_data.append({
                    'alpha': alpha_val, 'AR': ar_val, 
                    'CL': np.nan, 'CD': np.nan, 'L_D': np.nan
                })

    print("\nParameter sweep finished. Preparing plots...")

    # --- Plotting Results ---
    # For easier plotting, group results by AR
    grouped_results = {}
    for r in results_data:
        ar = r['AR']
        if ar not in grouped_results:
            grouped_results[ar] = {'alpha': [], 'CL': [], 'CD': [], 'L_D': []}
        grouped_results[ar]['alpha'].append(r['alpha'])
        grouped_results[ar]['CL'].append(r['CL'])
        grouped_results[ar]['CD'].append(r['CD'])
        grouped_results[ar]['L_D'].append(r['L_D'])

    # Plot 1: CL vs. Alpha for each AR
    plt.figure(figsize=(10, 6))
    for ar_val, data in grouped_results.items():
        plt.plot(data['alpha'], data['CL'], marker='o', linestyle='-', label=f'AR = {ar_val}')
    plt.title('CL vs. Angle of Attack (Alpha)')
    plt.xlabel('Alpha (degrees)')
    plt.ylabel('Lift Coefficient (CL)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Plot 2: L/D vs. Alpha for each AR
    plt.figure(figsize=(10, 6))
    for ar_val, data in grouped_results.items():
        plt.plot(data['alpha'], data['L_D'], marker='s', linestyle='-', label=f'AR = {ar_val}')
    plt.title('L/D Ratio vs. Angle of Attack (Alpha)')
    plt.xlabel('Alpha (degrees)')
    plt.ylabel('Lift-to-Drag Ratio (L/D)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Plot 3: Drag Polars (CL vs. CD) for each AR
    plt.figure(figsize=(10, 6))
    for ar_val, data in grouped_results.items():
        plt.plot(data['CD'], data['CL'], marker='x', linestyle='-', label=f'AR = {ar_val}')
    plt.title('Drag Polar (CL vs. CD)')
    plt.xlabel('Drag Coefficient (CD)')
    plt.ylabel('Lift Coefficient (CL)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    # Plot 4: Max L/D vs. AR
    max_ld_per_ar = []
    ars_for_plot4 = sorted(grouped_results.keys()) # Ensure ARs are plotted in order
    for ar_val in ars_for_plot4:
        data = grouped_results[ar_val]
        # Filter out nans before finding max, handle cases where all L/D are nan for an AR
        valid_lds = [ld for ld in data['L_D'] if not np.isnan(ld)]
        if valid_lds:
            max_ld_per_ar.append(np.max(valid_lds))
        else:
            max_ld_per_ar.append(np.nan) # Or 0, or skip this AR in plot

    plt.figure(figsize=(10, 6))
    plt.plot(ars_for_plot4, max_ld_per_ar, marker='D', linestyle='--')
    plt.title('Maximum L/D Ratio vs. Aspect Ratio (AR)')
    plt.xlabel('Aspect Ratio (AR)')
    plt.ylabel('Maximum L/D Ratio')
    plt.grid(True)
    plt.tight_layout()

    print("\nDisplaying plots...")
    plt.show()
    print("Script finished.")
```
