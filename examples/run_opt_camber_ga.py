"""Genetic algorithm optimization of wing camber distribution (assumed to be quadratic) 
to match target lift distribution.

Objective:
    Minimize sum of squared errors between actual and target local CL distribution

Design variables:
    Quadratic camber distribution: camber(y) = c + a*y^2
    - c: root camber (camber at y=0)
    - a: quadratic coefficient (constrained negative for downward-opening parabola)

    This ensures a continuous camber distribution across the span.

Flight condition:
    Velocity: 30 m/s
    Alpha: 1 degree

Approach:
    - Uses scipy.optimize.differential_evolution (genetic algorithm)
    - Uses NACA 4-digit airfoils specified directly in AVL geometry file
    - No separate .dat files needed - airfoils defined inline with NACA keyword
    - Evaluates fitness based on lift distribution match
"""

from optvl import OVLSolver
import numpy as np
from scipy.optimize import differential_evolution, NonlinearConstraint

# Configuration
NUM_SECTIONS  = 10
P_VALUE       = 0.4       # Location of max camber (40% chord) - must be 0.0-0.9 for NACA 4-digit
T_VALUE       = 0.09      # Thickness (9%) - must be 0.00-0.99 for NACA 4-digit
TEMP_AVL_FILE = "rectangular_wing_camber_temp.avl"
DEBUG         = False     # True for verbose output

# Trim conditions
TRIM_AOA      = 1.0       # Trim angle of attack (degrees)
TRIM_VELOCITY = 30.0      # Trim velocity (m/s)

# Wing twist (incidence angle distribution)
AOA_ROOT = 4.0      # Incidence angle at root (degrees)
AOA_TIP  = -1.0      # Incidence angle at tip (degrees) - negative is washout

# Section y-positions (half-span = 1.778m, 10 sections)
Y_TIP       = 1.778 # m, half span 
Y_POSITIONS = np.linspace(0, Y_TIP, NUM_SECTIONS)

# Quadratic camber distribution: camber(y) = c + a*y^2
# Design variable bounds
C_MIN, C_MAX = 0.0, 0.10      # Root camber: 0.5% to 5%
A_MIN, A_MAX = -0.04, 0.0     # Quadratic coefficient (negative for downward parabola)


def get_naca4_camber_slope_line(x_over_c, m, p):
    """
    Calculate NACA 4-digit camber line distribution.
    Calculate the slope (dyc/dx) of the NACA 4-digit camber line.

    Parameters:
        x_over_c: Array of x/c positions (0 to 1)
        m: Maximum camber as fraction of chord
        p: Location of maximum camber as fraction of chord

    Returns:
        Array of yc/c values (camber line ordinates)
        Array of dyc/dx values (camber line slopes)
    """
    yc = np.zeros_like(x_over_c)
    dyc_dx = np.zeros_like(x_over_c)

    for i, x in enumerate(x_over_c):
        if p == 0 or m == 0:
            yc[i] = 0.0
            dyc_dx[i] = 0.0
        elif x <= p:
            yc[i] = (m / p**2) * (2 * p * x - x**2)
            dyc_dx[i] = (2 * m / p**2) * (p - x)
        else:
            yc[i] = (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * x - x**2)
            dyc_dx[i] = (2 * m / (1 - p)**2) * (p - x)

    return yc, dyc_dx


def quadratic_camber_distribution(coeffs, y_positions):
    """
    Compute camber values at each span station. Assumes a quadratic distribution of camber.

    Parameters:
        coeffs: Array [c, a] where camber(y) = c + a*y^2
        y_positions: Array of span positions

    Returns:
        Array of camber values
    """
    c, a = coeffs
    cambers = c + a * y_positions**2
    return cambers


def linear_incidence_distribution(y_positions, aoa_root=AOA_ROOT, aoa_tip=AOA_TIP):
    """
    Compute linear incidence angle distribution from root to tip.

    Parameters:
        y_positions: Array of span positions
        aoa_root: Incidence angle at root (degrees)
        aoa_tip: Incidence angle at tip (degrees)

    Returns:
        Array of incidence angles (degrees)
    """
    # Linear interpolation: ainc(y) = aoa_root + (aoa_tip - aoa_root) * (y / y_tip)
    y_tip = y_positions[-1]
    incidences = aoa_root + (aoa_tip - aoa_root) * (y_positions / y_tip)
    return incidences


def create_avl_file_with_naca_airfoils(cambers, incidences, output_file):
    """
    Create an AVL file with NACA 4-digit airfoils specified directly.

    Parameters:
        cambers: Array of camber values (fraction of chord) for each section
        output_file: Path to output AVL file
    """


    avl_content = f"""AIRCRAFT
#======================================================
#------------------- Geometry File --------------------
#======================================================
# AVL Conventions
# SI Used: m, kg, etc
# Wing: 140 inches span, 9 inches chord, no taper
# Converted: 3.556m span, 0.2286m chord
# CAMBER OPTIMIZATION VARIANT - QUADRATIC DISTRIBUTION
# Twist: {AOA_ROOT:.1f}° root to {AOA_TIP:.1f}° tip

#Mach
0.0
#IYsym   IZsym   Zsym
 0       0       0
#Sref    Cref    b_wing
0.8127  0.2286  3.556
#Xref    Yref    Zref
0.0     0.0     0.0
# CDp
0.0

#======================================================
#--------------------- Main Wing ----------------------
#======================================================
SURFACE
Wing
#Nchordwise  Cspace   [Nspan   Sspace]
     8        1.0      20      -2.0
YDUPLICATE
0.0
SCALE
1.0  1.0  1.0
TRANSLATE
0.0  0.0  0.0
ANGLE
0.0

"""

    # Add sections with NACA airfoils and incidence angles
    for idx, (y_pos, camber, ainc) in enumerate(zip(Y_POSITIONS, cambers, incidences)):
        section_name = "Root" if idx == 0 else ("Tip" if idx == NUM_SECTIONS-1 else str(idx + 1))

        avl_content += f"""#------------------------------------------------------
# Section {idx + 1} ({section_name}) - Camber {camber*100:.2f}%, Ainc {ainc:.2f}°
SECTION
#Xle  Yle  Zle  |  Chord   Ainc   Nspan   Sspace
0.0    {y_pos:.4f}    0.0    0.2286   {ainc:.4f}    2      3

"""

    with open(output_file, 'w') as f:
        f.write(avl_content)


def get_target_lift_distribution(baseline_span, mean_chord):
    """
    Define the target lift distribution using Hunsaker's elliptical formula.
    Calculations must be in feet. 

    Reused from opt_rectangular_wing.py
    """
    span = baseline_span[-1] / 0.3048 * 2  # ft
    chord = mean_chord / 0.3048  # ft
    Sref = span * chord

    # lift distribution constant
    # =0 for elliptical lift distribution
    # =0.13564 for Hunsaker's elliptical lift distribution
    aconst = 0.13564 
    y_ft = baseline_span / 0.3048 # convert to ft
    eta = np.sqrt(1 - (2 * y_ft / span)**2)

    cl_dist = (4 * Sref) / (np.pi * span) * ((1 - 3 * aconst) * eta + 4 * aconst * eta**3)
    
    # trapz integral to get total lift
    total_lift = np.trapz(cl_dist, baseline_span/Y_TIP)

    return cl_dist, total_lift


def get_geometry(ovl):
    """Get geometric data from AVL strip analysis.

    Returns:
        baseline_span: Array of span positions (Y LE)
        mean_geometric_chord: Mean chord across all strips
    """
    strip_data    = ovl.get_strip_forces()
    wing_data     = strip_data["Wing"]
    baseline_span = np.array(wing_data['Y LE'])
    chord_values  = np.array(wing_data['chord'])
    mean_geometric_chord = np.mean(chord_values)

    return baseline_span, mean_geometric_chord


def get_current_cl_distribution(ovl):
    """Get the current CL distribution from AVL strip forces.

    Returns:
        cl_values: Array of local CL values at each strip
    """
    strip_data = ovl.get_strip_forces()
    wing_data = strip_data["Wing"]

    cl_values = np.array(wing_data['CL'])

    return cl_values


def evaluate_fitness(coeffs, ovl, target_cl):
    """
    Evaluate fitness using a reusable solver instance.

    Parameters:
        coeffs: Array [c, a] where camber(y) = c + a*y^2
        ovl: Pre-initialized OVLSolver instance
        target_cl: Target CL values

    Returns:
        float: Sum of squared errors (lower is better)
    """
    # Compute camber values from quadratic distribution
    cambers = quadratic_camber_distribution(coeffs, Y_POSITIONS)

    # Get the chordwise coordinate distribution from the solver
    # This is where AVL expects camber/thickness values
    x_dist = ovl.get_surface_param("Wing", "xasec")

    ovl.set_constraint('alpha', None, 1.0) # degrees
    ovl.set_parameter('velocity', 30.0) # m/sec

    # Get current casec, sasec, and ainsec arrays
    current_casec = ovl.get_surface_param("Wing", "casec")
    current_sasec = ovl.get_surface_param("Wing", "sasec")

    # Update camber and incidence for each section
    for i, m_val in enumerate(cambers):
        # Calculate NACA camber line for this section's x-distribution
        yc, dyc_dx = get_naca4_camber_slope_line(x_dist[i], m_val, P_VALUE)
        current_casec[i] = yc
        current_sasec[i] = dyc_dx

    # Push updated camber and incidence back to solver
    ovl.set_surface_param("Wing", "casec", current_casec)
    ovl.set_surface_param("Wing", "sasec", current_sasec)
    
    # Run analysis with updated geometry
    ovl.execute_run()
    current_cl = get_current_cl_distribution(ovl)

    # Calculate Sum of Squared Errors (SSE)
    errors = current_cl - target_cl
    sse = np.sum(errors**2)
        
    return float(sse)


# Global counter for progress tracking
eval_count = 0


def objective_wrapper(coeffs, ovl, target_cl):
    """Wrapper for objective function with progress printing."""
    global eval_count
    eval_count += 1

    sse = evaluate_fitness(coeffs, ovl, target_cl)

    if eval_count % 10 == 0:
        c, a = coeffs
        cambers = quadratic_camber_distribution(coeffs, Y_POSITIONS)
        print(f"  Eval {eval_count}: SSE = {sse:.6f}, "
              f"c={c*100:.3f}%, a={a*100:.4f}%/m^2, "
              f"root={cambers[0]*100:.2f}%, tip={cambers[-1]*100:.2f}%")

    return sse


def main():
    """Main optimization routine."""
    global eval_count

    print("=" * 60)
    print("CAMBER OPTIMIZATION USING GENETIC ALGORITHM")
    print("Quadratic Distribution: camber(y) = c + a*y^2")
    print("Using NACA 4-digit airfoils directly in AVL")
    print("=" * 60)

    # Run baseline to get span positions for target distribution
    print("\n=== BASELINE RUN ===")

    # Baseline: uniform 2% camber (c=0.02, a=0)
    baseline_coeffs  = np.array([0.02, 0.0])
    baseline_cambers = quadratic_camber_distribution(baseline_coeffs, Y_POSITIONS)
    incidences       = linear_incidence_distribution(Y_POSITIONS)
    
    geom_dict = {
        "title": "AIRCRAFT", # Aicraft name (MUST BE SET)
        "mach": 0.0, # Reference Mach number
        "iysym": 0, # y-symmetry settings
        "izsym": 0, # z-symmetry settings
        "zsym": 0.0, # z-symmetry plane
        "Sref": 0.8127, # Reference planform area
        "Cref": 0.2286, # Reference chord area
        "Bref": 3.556, # Reference span length
        "XYZref": np.array([0.0, 0, 0]), # Reference x,y,z position
        "CDp": 0.0, # Reference profile drag adjustment
        "surfaces": { # dictionary of surface dictionaries
            "Wing": {
                # General
                "num_sections": NUM_SECTIONS, # number of sections in surface
                "num_controls": np.array([0]*NUM_SECTIONS), # number of control surfaces assocaited with each section (in order)
                "num_design_vars": np.array([0]*NUM_SECTIONS), # number of AVL design variables assocaited with each section (in order)
                "component": 1,  # logical surface component index (for grouping interacting surfaces, see AVL manual)
                "yduplicate": 0.0,  # surface is duplicated over the ysymm plane
                "xles": np.zeros(NUM_SECTIONS),  # leading edge cordinate vector(x component)
                "yles": np.linspace(0, Y_TIP, NUM_SECTIONS),  # leading edge cordinate vector(y component)
                "zles": np.zeros(NUM_SECTIONS),  # leading edge cordinate vector(z component)
                "chords": np.ones(NUM_SECTIONS)*0.2286,  # chord length vector
                "aincs": incidences,  # incidence angle vector
                # NACA
                'naca' : np.array(['2412']*NUM_SECTIONS), # 4-digit NACA airfoil
                # Paneling
                "nchordwise": 8,  # number of chordwise horseshoe vortice s placed on the surface
                "cspace": 1.0,  # chordwise vortex spacing parameter
                "nspan": 20,  # number of spanwise horseshoe vortices placed on the entire surface
                "sspace": -2.0,  # spanwise vortex spacing parameter for entire surface
                "nspans": np.array([5]*NUM_SECTIONS),  # number of spanwise elements vector
                "sspaces": np.array([3.0]*NUM_SECTIONS),  # spanwise spacing vector (for each section)
                "use surface spacing": 1,  # surface spacing set for the entire surface (known as LSURFSPACING in AVL)
        },
        }
    }

    # Initialize solver ONCE - this will be reused throughout optimization
    ovl = OVLSolver(input_dict=geom_dict, debug=DEBUG)
    ovl.set_constraint('alpha', None, TRIM_AOA)
    ovl.set_parameter('velocity', TRIM_VELOCITY)
    ovl.execute_run()
    
    baseline_forces = ovl.get_total_forces()
    baseline_span, mean_chord = get_geometry(ovl)
    baseline_cl = get_current_cl_distribution(ovl)

    print(f"Baseline CL (total): {baseline_forces['CL']:.6f}")
    print(f"Baseline CD (total): {baseline_forces['CD']:.6f}")
    print(f"Baseline L/D: {baseline_forces['CL']/baseline_forces['CD']:.6f}")

    # Get target distribution
    target_cl, total_lift = get_target_lift_distribution(baseline_span, mean_chord)

    print("\n=== Target Lift Distribution ===")
    print("Span (m)  | Target CL")
    print("-" * 25)
    for span, cl in zip(baseline_span, target_cl):
        print(f"{span:8.4f}  | {cl:8.4f}")

    baseline_sse = evaluate_fitness(baseline_coeffs, ovl, target_cl)
    print(f"\nBaseline SSE (2% uniform camber): {baseline_sse:.6f}")

    # Set up bounds for quadratic coefficients [c, a]
    # a is constrained to be negative (parabola opens downward)
    bounds = [
        (C_MIN, C_MAX),  # c: root camber
        (A_MIN, A_MAX),  # a: quadratic coefficient (negative)
    ]

    # Define constraint: c + a*Y_TIP^2 >= 0 (ensures positive camber at tip)
    def tip_camber_constraint(x):
        """Ensure camber at tip is non-negative: c + a*Y_TIP^2 >= 0"""
        c, a = x
        return c + a * Y_TIP**2

    constraint = NonlinearConstraint(tip_camber_constraint, 1e-4, np.inf)  # Minimum 0.1% camber at tip

    print("\n=== STARTING GENETIC ALGORITHM OPTIMIZATION ===")
    print("Design variables: 2 quadratic coefficients [c, a]")
    print(f"  c (root camber):      {C_MIN*100:.1f}% to {C_MAX*100:.1f}%")
    print(f"  a (quadratic coeff):  {A_MIN*100:.2f}%/m^2 to {A_MAX*100:.2f}%/m^2")
    print(f"Constraint: c + a*Y_TIP^2 >= 0.001 (ensures positive camber at tip)")
    print()

    eval_count = 0

    # Run differential evolution (genetic algorithm)
    result = differential_evolution(
        func=lambda x: objective_wrapper(x, ovl, target_cl),
        bounds=bounds,
        constraints=[constraint],
        strategy='best1bin',
        maxiter=50,
        popsize=15,
        tol=1e-7,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        disp=True,
        workers=1,  # Use 1 worker since we're reusing solver instance
        updating='immediate'
    )

    print("\n" + "=" * 60)
    print("OPTIMIZATION RESULTS")
    print("=" * 60)

    optimal_coeffs = result.x
    c_opt, a_opt = optimal_coeffs
    optimal_cambers = quadratic_camber_distribution(optimal_coeffs, Y_POSITIONS)
    optimal_cambers_pct = optimal_cambers * 100

    print(f"\nOptimization converged: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of evaluations: {result.nfev}")
    print(f"Final SSE: {result.fun:.6f}")
    print(f"Improvement: {((baseline_sse - result.fun) / baseline_sse * 100):.2f}%")

    print("\nOptimal Quadratic Coefficients:")
    print(f"  c (root camber) = {c_opt*100:.4f}%")
    print(f"  a (quadratic)   = {a_opt*100:.4f}%/m^2")
    print(f"\nCamber distribution: camber(y) = {c_opt*100:.4f} + {a_opt*100:.4f}*y^2  [%]")

    print("\nResulting Camber at Each Section:")
    print("Section  | Span (m) | Camber (%)")
    print("-" * 40)
    for i, (y, camber) in enumerate(zip(Y_POSITIONS, optimal_cambers_pct)):
        print(f"   {i+1:2d}    | {y:8.4f} | {camber:8.3f}")

    # Generate final configuration and analyze
    print("\n=== FINAL CONFIGURATION ANALYSIS ===")

    # Save final configuration file for reference
    create_avl_file_with_naca_airfoils(optimal_cambers, incidences, "rectangular_wing_camber_optimized.avl")

    # Reuse existing solver by updating its geometry
    x_dist = ovl.get_surface_param("Wing", "xasec")
    current_casec = ovl.get_surface_param("Wing", "casec")
    current_sasec = ovl.get_surface_param("Wing", "sasec")
    
    for i, m_val in enumerate(optimal_cambers):
        yc, dyc_dx = get_naca4_camber_slope_line(x_dist[i], m_val, P_VALUE)
        current_casec[i] = yc
        current_sasec[i] = dyc_dx

    ovl.set_surface_param("Wing", "casec", current_casec)
    ovl.set_surface_param("Wing", "sasec", current_sasec)

    ovl.execute_run()

    final_cl = get_current_cl_distribution(ovl)
    final_forces = ovl.get_total_forces()

    print(f"\nFinal CL (total): {final_forces['CL']:.6f}")
    print(f"Final CD (total): {final_forces['CD']:.6f}")
    print(f"Final L/D: {final_forces['CL']/final_forces['CD']:.6f}")

    # Compare distributions
    print("\nLift Distribution Comparison:")
    print("Span (m)  | Target CL | Final CL | Error")
    print("-" * 50)

    for span, target, final in zip(baseline_span, target_cl, final_cl):
        error = final - target
        print(f"{span:8.4f}  | {target:9.4f} | {final:8.4f} | {error:+7.4f}")

    print(f"\nOptimized AVL file saved to: rectangular_wing_camber_optimized.avl")
    print("NACA airfoils specified directly in AVL file (no separate .dat files needed)")

    # ######## PLOT RESULTS ########
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Plot 1: Lift distribution comparison
    ax1 = axes[0]
    ax1.plot(baseline_span, target_cl, 'o-', label='Target', linewidth=2, markersize=8)
    ax1.plot(baseline_span, baseline_cl, 's--', label='Baseline (2% uniform)', linewidth=1.5, markersize=6)
    ax1.plot(baseline_span, final_cl, '^-', label='Optimized', linewidth=1.5, markersize=6)
    ax1.set_xlabel('Span Position (m)', fontsize=12)
    ax1.set_ylabel('Local CL', fontsize=12)
    ax1.set_title('Lift Distribution', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Continuous quadratic camber distribution
    ax2 = axes[1]
    # Plot continuous curve
    y_continuous = np.linspace(0, Y_TIP, 100)
    camber_continuous = quadratic_camber_distribution(optimal_coeffs, y_continuous) * 100
    ax2.plot(y_continuous, camber_continuous, 'b-', linewidth=2, label='Quadratic fit')

    # Plot discrete section values
    ax2.scatter(Y_POSITIONS, optimal_cambers_pct, s=80, c='steelblue',
               edgecolors='navy', zorder=5, label='Section values')
    ax2.axhline(y=2.0, color='red', linestyle='--', label='Baseline (2%)')
    ax2.set_xlabel('Span Position (m)', fontsize=12)
    ax2.set_ylabel('Camber (%)', fontsize=12)
    ax2.set_title(f'Quadratic Camber: {c_opt*100:.2f} + {a_opt*100:.2f}y$^2$',
                 fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Camber lines at key sections
    ax3 = axes[2]
    sections_to_plot = [0, 4, 9]  # Root, mid, tip
    colors = ['blue', 'green', 'red']
    labels = ['Root', 'Mid', 'Tip']
    x_camber = np.linspace(0, 1, 100)  # x/c positions
    for sec_idx, color, label in zip(sections_to_plot, colors, labels):
        yc, _ = get_naca4_camber_slope_line(x_camber, optimal_cambers[sec_idx], P_VALUE)
        # Offset vertically for visibility
        offset = (sections_to_plot.index(sec_idx)) * 0.03
        ax3.plot(x_camber, yc + offset, color=color, linewidth=2,
                label=f'{label} ({optimal_cambers_pct[sec_idx]:.2f}%)')
    ax3.set_xlabel('x/c', fontsize=12)
    ax3.set_ylabel('Camber (y/c, offset)', fontsize=12)
    ax3.set_title('Optimized Camber Lines', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('opt_camber_ga_results.png', dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: opt_camber_ga_results.png")
    plt.show()

if __name__ == "__main__":
    main()
