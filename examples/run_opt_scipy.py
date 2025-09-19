"""A scipy based optimization to trim an aircraft for an aicraft using optvl"""
import numpy as np
from scipy.optimize import minimize
from optvl import OVLSolver

ovl_solver = OVLSolver(geo_file="aircraft.avl", debug=False)

# setup OptVL 
ovl_solver.set_parameter("Mach", 0.0)

# obj-start
# Define your custom objective function with outputs from OptVL
def objective_function(x):
    ovl_solver.set_control_deflection("Elevator", x[0])
    ovl_solver.set_surface_params({"Wing":{"aincs":x[1:]}})
    
    ovl_solver.execute_run()
    cd = ovl_solver.get_total_forces()['CD']
    print(x, cd)

    return cd
# obj-end
# objgrad-start
def objective_gradient(x):
    # Partial derivatives of the objective_function
    
    # we are trusting that the design variables have already been applied 
    # and propogated through by the objective_function. 
    
    ovl_solver.set_control_deflection("Elevator", x[0])
    ovl_solver.set_surface_params({"Wing":{"aincs":x[1:]}})
    
    ovl_solver.execute_run()
    
    sens = ovl_solver.execute_run_sensitivities(['CD'])
    dcd_dele = sens['CD']['Elevator']
    dcd_daincs = sens['CD']['Wing']['aincs']

    # concatinate the two and return the derivs
    return np.concatenate(([dcd_dele], dcd_daincs))
# objgrad-end

# eq-start
# Define equality constraint: h(x) = 0
def eq_constraint(x):
    ovl_solver.set_control_deflection("Elevator", x[0])
    ovl_solver.set_surface_params({"Wing":{"aincs":x[1:]}})

    ovl_solver.execute_run()

    # the objective must always be run first
    coeff = ovl_solver.get_total_forces()
    
    cl_target = 1.5
    cm_target = 0.0

    cl_con = coeff['CL'] - cl_target
    cm_con = coeff['Cm'] - cm_target

    return np.array([cl_con, cm_con])
# eq-end
# eqgrad-start
# Define the gradient of the equality constraint
def eq_constraint_jac(x):
    sens = ovl_solver.execute_run_sensitivities(['CL', 'Cm'])
    dcl_dele = sens['CL']['Elevator']
    dcl_daincs = sens['CL']['Wing']['aincs']
    dcl_dx = np.concatenate(([dcl_dele], dcl_daincs))
    
    dcm_dele = sens['Cm']['Elevator']
    dcm_daincs = sens['Cm']['Wing']['aincs']
    dcm_dx = np.concatenate(([dcm_dele], dcm_daincs))

    # concatinate the two and return the derivs
    return np.array([dcl_dx, dcm_dx])
# eqgrad-end

num_sec = 5
# Initial guess for the variables
x0 = np.zeros(1+num_sec)

# Define constraints with their gradients
constraints = [
    {'type': 'eq', 'fun': eq_constraint, 'jac': eq_constraint_jac},   # Equality constraint
]

# Call the minimize function with constraints and their gradients
result = minimize(
    fun=objective_function,       # The objective function to minimize
    jac=objective_gradient,       # The gradient function of the objective
    x0=x0,                     # Initial guess
    constraints=constraints,   # Constraints with gradients
    method='SLSQP',            # Optimization method that supports constraints
    options={'disp': True},     # Display convergence messages
    tol=1e-10
)

# Print the result
print("Optimization result:")
opt_x = result.x
print(f"Elevator: {opt_x[0]}")
print(f"Wing twist: {opt_x[1:]}")
print(f"Final CD: {result.fun}")          # Final value of the objective function
