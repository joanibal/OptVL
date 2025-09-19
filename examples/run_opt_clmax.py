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
from optvl import OVLSolver
import matplotlib.pyplot as plt

# Initialize OptVL solver
geo_file = "aircraft.avl" 
ovl_solver = OVLSolver(geo_file=geo_file, debug=False)


# Set Mach number for the analysis
ovl_solver.set_parameter("Mach", 0.0)

def objective_function(dvs):
    """
    Calculates the negative of the total lift coefficient (CL) for a given angle of attack.
    Scipy's minimize function tries to minimize this value.
    """
    alpha = dvs[0] 
    ovl_solver.set_variable("alpha", alpha) 
    def_elev = dvs[1]
    ovl_solver.set_control_deflection("Elevator", def_elev) 
    
    ovl_solver.execute_run()
    
    total_forces = ovl_solver.get_total_forces()
    cl = total_forces['CL']
    print(f"Objective: Alpha={alpha:.4f} deg, CL={cl:.4f}")
    
    return -cl # Scipy minimizes, so return negative CL to maximize CL


def cl_con(dvs):
    """
    Ensures that no sectional lift coefficient (Cl) exceeds 1.75.
    Returns a value >= 0 if the constraint is met (1.75 - max_cl_sectional).
    """
    alpha = dvs[0]
    ovl_solver.set_variable("alpha", alpha)
    def_elev = dvs[1]
    ovl_solver.set_control_deflection("Elevator", def_elev) 
    ovl_solver.execute_run()
    
    strip_forces = ovl_solver.get_strip_forces()
    sectional_cls = strip_forces['Wing']['CL']
    
    sectional_cls_limit = np.ones_like(sectional_cls)*1.75
    
    con_values = sectional_cls_limit - sectional_cls
    
    return con_values

def cm_con(dvs):
    alpha = dvs[0]
    ovl_solver.set_variable("alpha", alpha)
    def_elev = dvs[1]
    ovl_solver.set_control_deflection("Elevator", def_elev) 
    ovl_solver.execute_run()
    
    cm = ovl_solver.get_total_forces()['Cm']
    
    return cm

# Initial guess for design variables
x0 = np.array([10.0, 0]) 

constraints = [{'type': 'ineq', 'fun': cl_con}, {'type': 'eq', 'fun': cm_con}]

# Run the optimization using SLSQP
# - If we don't supply gradients then the optimizer will use finite difference
# - Since we only have a few design variables it is fine to use finite difference in this case.
print("\nStarting optimization...")
result = minimize(
    fun=objective_function,    # Objective function to minimize
    x0=x0,                     # Initial guess
    method='SLSQP',            # Sequential Least Squares Programming
    constraints=constraints,   # Constraints definition
    options={'disp': True, 'ftol': 1e-6} # Display convergence and set tolerance
)

if result.success:
    print("Optimization successful!")
else:
    print("Something went wrong in the optimization!")

# --- Results and Verification ---

optimized_alpha = result.x[0]
optimized_def_elev = result.x[1]
calculated_cl_max = -result.fun 

print(f"Optimized Angle of Attack: {optimized_alpha:.4f} degrees")
print(f"Optimized Elevator deflection: {optimized_def_elev:.4f} degrees")

print(f"Calculated CL_max from optimizer: {calculated_cl_max:.4f}") # Objective was -CL

# Verification run at the optimized alpha to double-check AVL results
print("\nVerification run at optimized alpha:")
ovl_solver.set_variable("alpha", optimized_alpha)
ovl_solver.set_control_deflection("Elevator", optimized_def_elev)


# Recalculate max sectional Cl at the optimal alpha for verification
ovl_solver.execute_run()
total_forces_verify = ovl_solver.get_total_forces()
cl_verify = total_forces_verify['CL']
cd_verify = total_forces_verify['CD'] 
cm_verify = total_forces_verify['Cm']
strip_forces = ovl_solver.get_strip_forces()
sectional_cls = strip_forces['Wing']['CL']
y_le = strip_forces['Wing']['Y LE']


print(f"  Total CL from verification run: {cl_verify:.4f}")
print(f"  Total CD from verification run: {cd_verify:.4f}")
print(f"  Total Cm from verification run: {cm_verify:.4f}")
print(f"  Maximum sectional Cl from verification run: {np.max(sectional_cls):.4f} (Constraint: <= 1.75)")

# --- Plotting Span-wise Cl Distribution ---
plt.figure(figsize=(10, 6))
        
plt.plot(y_le, sectional_cls, label='sectional Cl')
plt.axhline(y=1.75, color='r', linestyle='--', label='Max Sectional Cl (1.75)')

plt.xlabel("Span-wise location (Y coordinate of LE strip)")
plt.ylabel("Sectional Lift Coefficient (Cl)")
plt.title(f"Span-wise Cl Distribution at Alpha = {optimized_alpha:.2f} deg (CL_max = {cl_verify:.3f})")
plt.legend()
plt.grid(True)
plt.show()
