"""
Example script to find the maximum lift coefficient (CL_max) of an aircraft
using the critical section method with Scipy's SLSQP optimizer and OptVL.

In the critical section method, the aircraft is considered to begin stalling,
and its current CL is taken as CL_max, when any chord-wise section's local
lift coefficient (Cl_sectional) exceeds a predefined maximum value (Cl_max_sectional).
This script varies the aircraft's angle of attack (alpha) to maximize the
total CL while ensuring no sectional Cl goes above 1.75.
"""
import numpy as np
from scipy.optimize import minimize
from optvl import OVLSolver # Adjusted import based on typical structure and OVLSolver definition
import matplotlib.pyplot as plt

# Initialize OptVL solver
# Assumes 'aircraft.avl' is in the 'examples/' directory and the script is run from the repo root.
geo_file = "examples/aircraft.avl" 
ovl_solver = OVLSolver(geo_file=geo_file, debug=False)


# Set Mach number for the analysis
ovl_solver.set_parameter("Mach", 0.0)

print("OVLSolver initialized and Mach number set.")


def objective_function(alpha_deg):
    """
    Calculates the negative of the total lift coefficient (CL) for a given angle of attack.
    Scipy's minimize function tries to minimize this value.
    """
    alpha = float(alpha_deg) # Ensure alpha is float
    # print(f"Objective: Setting alpha to {alpha:.4f} deg")
    ovl_solver.set_constraint("alpha", alpha) # Alpha is the primary variable for AVL
    
    # Ensure other controls (elevator, etc.) are fixed or at default for CL_max analysis.
    # This example assumes default control surface deflections as per the .avl file.
    
    ovl_solver.execute_run()
    
    total_forces = ovl_solver.get_total_forces()
    cl = total_forces['CL']
    # print(f"Objective: Alpha={alpha:.4f} deg, CL={cl:.4f}")
    
    return -cl # Scipy minimizes, so return negative CL to maximize CL


def constraint_function(alpha_deg):
    """
    Ensures that no sectional lift coefficient (Cl) exceeds 1.75.
    Returns a value >= 0 if the constraint is met (1.75 - max_cl_sectional).
    """
    alpha = float(alpha_deg) # Ensure alpha is a float
    # print(f"Constraint: Setting alpha to {alpha:.4f} deg")
    ovl_solver.set_constraint("alpha", alpha)
    ovl_solver.execute_run()
    
    strip_forces = ovl_solver.get_strip_forces()
    all_sectional_cls = []
    
    if not strip_forces:
        # Handle missing strip_forces (e.g., run error, no lifting surfaces)
        # Return large negative for severe constraint violation.
        print("Warning: No strip forces found. Check AVL execution and geometry.")
        return -100.0 
        
    for surface_name in strip_forces:
        if 'CL strip' in strip_forces[surface_name] and strip_forces[surface_name]['CL strip'].size > 0:
            all_sectional_cls.append(strip_forces[surface_name]['CL strip'])
        # else:
            # print(f"Warning: 'CL strip' not found or empty for surface {surface_name}")

    if not all_sectional_cls:
        # Handle no sectional Cls (no lifting surfaces with 'CL strip' data).
        # Return large negative for severe constraint violation.
        print("Warning: No sectional Cls collected from any surface.")
        return -100.0

    # Filter out potential NaNs or Infs if necessary, though AVL should provide valid numbers.
    # The line below is an example of such filtering if it were needed.
    # concatenated_cls = np.concatenate([cls[np.isfinite(cls)] for cls in all_sectional_cls if cls.size > 0])
    
    # Ensure all collected sectional Cl arrays are not empty before concatenation
    valid_cls_arrays = [cls for cls in all_sectional_cls if cls.size > 0]
    
    if not valid_cls_arrays:
        print("Warning: All sectional Cl arrays were empty after initial check.")
        return -100.0 # Treat as violation, implies no valid data to process

    concatenated_cls = np.concatenate(valid_cls_arrays)

    if concatenated_cls.size == 0:
        # This should ideally not be reached if valid_cls_arrays check is robust
        print("Warning: Concatenated sectional Cls array is empty.")
        return -100.0 # Treat as violation

    max_cl_sectional = np.max(concatenated_cls)
    # print(f"Constraint: Alpha={alpha:.4f} deg, Max Sectional Cl={max_cl_sectional:.4f}")
    
    # Constraint: 1.75 - max_cl_sectional >= 0
    return 1.75 - max_cl_sectional


def objective_gradient(alpha_deg):
    """
    Calculates the gradient of the objective function (-CL) with respect to alpha.
    This is - (dCL/dalpha).
    """
    alpha = float(alpha_deg) # Ensure alpha is float
    # print(f"Gradient: Setting alpha to {alpha:.4f} deg")
    ovl_solver.set_constraint("alpha", alpha)
    
    # Execute run to update AVL state before calculating sensitivities
    ovl_solver.execute_run() 
    
    sens = ovl_solver.execute_run_sensitivies(funcs=['CL']) # Request sensitivity of CL
    
    dcl_dalpha = 0.0 # Default if sensitivity not found
    if 'CL' in sens and 'alpha' in sens['CL']:
        dcl_dalpha = sens['CL']['alpha']
    else:
        print(f"Warning: Sensitivity dCL/dalpha not found for alpha={alpha:.4f}. Returning 0.0 gradient.")
        print(f"Sensitivity dict: {sens}")

    # print(f"Gradient: Alpha={alpha:.4f} deg, dCL/dAlpha={dcl_dalpha:.4f}")
    
    return -dcl_dalpha # Gradient of -CL is -(dCL/dalpha)

# Initial guess for alpha (degrees)
x0 = np.array([2.0]) 

# Bounds for alpha (degrees); adjust as needed for the specific aircraft
alpha_bounds = (-5.0, 20.0)

# Define constraints for the optimizer
# The constraint is inequality: constraint_function(alpha) >= 0
constraints = [{'type': 'ineq', 'fun': constraint_function}]

# Run the optimization using SLSQP
print("\nStarting optimization...")
result = minimize(
    fun=objective_function,    # Objective function to minimize
    x0=x0,                     # Initial guess
    jac=objective_gradient,    # Jacobian (gradient) of the objective function
    method='SLSQP',            # Sequential Least Squares Programming
    bounds=[alpha_bounds],     # Bounds for the optimization variable (alpha)
    constraints=constraints,   # Constraints definition
    options={'disp': True, 'ftol': 1e-6} # Display convergence and set tolerance
)

# --- Results and Verification ---
if result.success:
    print("\nOptimization successful!")
    optimized_alpha = result.x[0]
    calculated_cl_max = -result.fun 

    print(f"Optimized Angle of Attack (alpha_opt): {optimized_alpha:.4f} degrees")
    print(f"Calculated CL_max from optimizer: {calculated_cl_max:.4f}") # Objective was -CL

    # Verification run at the optimized alpha to double-check AVL results
    print("\nVerification run at optimized alpha:")
    ovl_solver.set_constraint("alpha", float(optimized_alpha))
    ovl_solver.execute_run()
    
    total_forces_verify = ovl_solver.get_total_forces()
    cl_verify = total_forces_verify['CL']
    cd_verify = total_forces_verify['CD'] 
    cm_verify = total_forces_verify['CM']

    # Recalculate max sectional Cl at the optimal alpha for verification
    # constraint_function itself calls execute_run, so AVL is run again here.
    constraint_val_at_opt = constraint_function(optimized_alpha) 
    max_sectional_cl_at_opt = 1.75 - constraint_val_at_opt
    
    print(f"  Total CL from verification run: {cl_verify:.4f}")
    print(f"  Total CD from verification run: {cd_verify:.4f}")
    print(f"  Total CM from verification run: {cm_verify:.4f}")
    print(f"  Maximum sectional Cl from verification run: {max_sectional_cl_at_opt:.4f} (Constraint: <= 1.75)")

    # Check consistency between optimizer result and verification run
    if not np.isclose(calculated_cl_max, cl_verify, atol=1e-4):
        print(f"  Warning: Optimizer CL_max ({calculated_cl_max:.4f}) and verification CL ({cl_verify:.4f}) differ.")
    
    # Check if constraint is met at the found optimum (with a small tolerance)
    if max_sectional_cl_at_opt > 1.75 + 1e-4: 
        print(f"  Warning: Constraint VIOLATED at optimal alpha. Max sectional Cl = {max_sectional_cl_at_opt:.4f}")

    # --- Plotting Span-wise Cl Distribution ---
    print("\nGenerating span-wise Cl distribution plot...")
    
    # Get strip forces again for the optimized alpha if not already available
    # Note: The verification run already set alpha and ran the solver.
    # We can call get_strip_forces directly.
    strip_forces_plot = ovl_solver.get_strip_forces()
    
    plt.figure(figsize=(10, 6))
    
    # Plot for specific surfaces, e.g., "Wing"
    # Adjust surface_names_to_plot if your AVL file uses different names or you want more surfaces
    surface_names_to_plot = ["Wing"] 
    
    plotted_something = False
    for surface_name in surface_names_to_plot:
        if surface_name in strip_forces_plot and 'Y LE' in strip_forces_plot[surface_name] and 'CL strip' in strip_forces_plot[surface_name]:
            y_le = strip_forces_plot[surface_name]['Y LE']
            cl_strip = strip_forces_plot[surface_name]['CL strip']
            
            if y_le.size > 0 and cl_strip.size > 0 and y_le.size == cl_strip.size:
                # Sort by Y LE for a clean plot, especially if data isn't ordered
                sorted_indices = np.argsort(y_le)
                plt.plot(y_le[sorted_indices], cl_strip[sorted_indices], marker='o', linestyle='-', label=f'{surface_name} Cl distribution')
                plotted_something = True
            else:
                print(f"  Skipping plot for {surface_name}: Y LE or CL strip data is empty or mismatched.")
        else:
            print(f"  Skipping plot for {surface_name}: Data not found in strip_forces_plot.")

    if plotted_something:
        plt.axhline(y=1.75, color='r', linestyle='--', label='Max Sectional Cl (1.75)')
        plt.xlabel("Span-wise location (Y coordinate of LE strip)")
        plt.ylabel("Sectional Lift Coefficient (Cl)")
        plt.title(f"Span-wise Cl Distribution at Alpha = {optimized_alpha:.2f} deg (CL_max = {cl_verify:.3f})")
        plt.legend()
        plt.grid(True)
        plt.show()
    else:
        print("  Could not plot span-wise Cl distribution: No suitable surface data found or plotted.")

else:
    print("\nOptimization failed.")
    print(f"  Message: {result.message}")
    print(f"  Number of iterations: {result.nit}")
    print(f"  Function evaluations: {result.nfev}")
    if hasattr(result, 'njev'): # Jacobian evaluations
        print(f"  Jacobian evaluations: {result.njev}")
    
    # Display state at the last attempted alpha if optimization failed
    last_alpha_attempt = result.x[0]
    print(f"\nState at last attempted alpha = {last_alpha_attempt:.4f} degrees:")
    
    obj_at_last = objective_function(last_alpha_attempt) # Runs solver
    print(f"  Objective function value (-CL): {obj_at_last:.4f} (CL: {-obj_at_last:.4f})")
    
    constraint_val_at_last = constraint_function(last_alpha_attempt) # Runs solver
    max_cl_sect_at_last = 1.75 - constraint_val_at_last
    print(f"  Constraint function value (1.75 - max_sectional_cl): {constraint_val_at_last:.4f}")
    print(f"  Max sectional Cl: {max_cl_sect_at_last:.4f}")

print("\nEnd of script.")
