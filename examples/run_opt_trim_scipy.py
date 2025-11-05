"""A scipy based optimization to trim an aircraft for an aicraft using optvl"""

import numpy as np
from scipy.optimize import minimize
from optvl import OVLSolver

ovl_solver = OVLSolver(geo_file="../geom_files/aircraft.avl", debug=False)

# setup OptVL
ovl_solver.set_parameter("Mach", 0.0)
# ovl_solver.set_variable("alpha", 5.0)


# Define your custom objective function with outputs from OptVL
def custom_function(x):
    ovl_solver.set_control_deflection("Elevator", x[0])
    ovl_solver.set_variable("alpha", x[1])

    ovl_solver.execute_run()
    cd = ovl_solver.get_total_forces()["CD"]
    return cd


# Define the gradient (Jacobian) of the objective function
def custom_gradient(x):
    # Partial derivatives of the custom_function
    sens = ovl_solver.execute_run_sensitivities(["CD"])
    dcd_dele = sens["CD"]["Elevator"]
    dcd_dalpha = sens["CD"]["alpha"]
    # concatinate the two and return the derivs
    return np.array([dcd_dele, dcd_dalpha])


cl_target = 1.5
cm_target = 0.0


# Define equality constraint: h(x) = 0
def eq_constraint(x):
    ovl_solver.set_control_deflection("Elevator", x[0])
    ovl_solver.set_variable("alpha", x[1])
    ovl_solver.execute_run()

    # the objective must always be run first
    coeff = ovl_solver.get_total_forces()

    cl_con = coeff["CL"] - cl_target
    cm_con = coeff["Cm"] - cm_target
    return np.array([cl_con, cm_con])


# Define the gradient of the equality constraint
def eq_constraint_jac(x):
    sens = ovl_solver.execute_run_sensitivities(["CL", "Cm"])
    dcl_dele = sens["CL"]["Elevator"]
    dcl_dalpha = sens["CL"]["alpha"]

    dcm_dele = sens["Cm"]["Elevator"]
    dcm_dalpha = sens["Cm"]["alpha"]
    # concatinate the two and return the derivs
    return np.array([[dcl_dele, dcl_dalpha], [dcm_dele, dcm_dalpha]])


# Initial guess for the variables
x0 = np.array([0.0, 0.0])

# Define constraints with their gradients
constraints = [
    {"type": "eq", "fun": eq_constraint, "jac": eq_constraint_jac},  # Equality constraint
]

# Call the minimize function with constraints and their gradients
result = minimize(
    fun=custom_function,  # The objective function to minimize
    jac=custom_gradient,  # The gradient function of the objective
    x0=x0,  # Initial guess
    constraints=constraints,  # Constraints with gradients
    method="SLSQP",  # Optimization method that supports constraints
    options={"disp": True},  # Display convergence messages
)

# Print the result
print("Optimization result:")
opt_x = result.x
print(f"Elevator:{opt_x[0]}, AoA:{opt_x[1]}")
print(f"Final CD:{result.fun}")  # Final value of the objective function
