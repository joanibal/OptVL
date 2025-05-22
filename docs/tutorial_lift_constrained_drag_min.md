# Tutorial: Lift-Constrained Drag Minimization of a Wing

This tutorial demonstrates how to perform a common aerodynamic optimization task: minimizing the drag coefficient (CD) of a wing while maintaining a specific lift coefficient (CL). This is a classic problem in aircraft design, often relevant for optimizing cruise performance.

## 1. Problem Definition

*   **Objective:** Minimize the drag coefficient (CD) of the wing.
*   **Constraint:** Maintain a fixed lift coefficient (CL) of 0.5. This value is representative of typical cruise conditions for many aircraft.
*   **Aircraft:** We will use the simple rectangular wing geometry defined in `examples/rectangle.avl`. This provides a straightforward starting point for understanding the optimization process.
*   **Flight conditions:** The analysis will be performed at fixed flight conditions relevant to cruise (e.g., Mach 0.75, altitude of 10,000 meters). Specific values will be set in the OptVL analysis.

## 2. Design Variables

For this tutorial, we will use the following design variables that OptVL can modify:

*   **Angle of Attack (alpha):**
    *   This is the primary variable the optimizer will use to meet the target lift coefficient.
    *   Range: Typically, -5 to 10 degrees. The optimizer will find the specific alpha needed.
    *   Initial Value: Around 2.0 degrees.

*   **Wing Aspect Ratio (AR):**
    *   Aspect ratio (AR = span^2 / area) is a key parameter influencing induced drag. We will assume OptVL allows modification of span while keeping the wing area constant (or vice-versa, by adjusting chord). If direct AR modification isn't straightforward, one might choose span or chord as a variable, with area being recalculated or constrained. For this conceptual tutorial, we'll assume AR can be varied.
    *   Range: e.g., 6 to 12.
    *   Initial Value: The AR defined in `examples/rectangle.avl` (e.g., if it's 8, we can start there).

*   **Linear Wing Twist (Washout):**
    *   This refers to a linear variation in the angle of incidence of the wing sections from the root to the tip. Negative twist (washout) means the tip is at a lower angle of incidence than the root, which can improve spanwise lift distribution and reduce induced drag.
    *   We'll define this as the twist angle at the tip (e.g., root twist is 0, tip twist is the design variable).
    *   Range: e.g., -5 degrees to 0 degrees.
    *   Initial Value: 0 degrees (no twist).

## 3. Optimization Setup Overview

We will use an optimization algorithm (this tutorial will later show an example using Scipy) to adjust the design variables: angle of attack, aspect ratio, and wing twist.

The optimizer's task is to explore different combinations of these variables. For each combination, OptVL will calculate the resulting CL and CD. The optimizer then aims to:
1.  Find the set of design variables that produces the lowest possible CD.
2.  While simultaneously ensuring that the CL for that design is exactly 0.5 (or very close to it).

OptVL plays a crucial role by providing not just the CL and CD values, but also the derivatives (sensitivities) of CL and CD with respect to each design variable. This derivative information helps the optimizer make intelligent decisions about how to adjust the variables to move towards the optimal solution more efficiently.

## 4. OptVL's Role

In this optimization loop, OptVL will be responsible for:

1.  **Loading Geometry:** Reading the baseline aircraft geometry from the `examples/rectangle.avl` file.
2.  **Setting Flight Conditions:** Establishing the specified Mach number and altitude for the analysis.
3.  **Applying Design Variables:**
    *   Setting the current angle of attack (alpha) for the analysis.
    *   Modifying the wing geometry based on the current aspect ratio. This might involve scaling the span and/or chord appropriately while respecting any constraints (e.g., constant area if specified).
    *   Applying the specified linear twist to the wing sections.
4.  **Executing Analysis:** Running the AVL core routines to perform the aerodynamic calculations for the modified geometry and flight state.
5.  **Returning Results:** Providing the calculated CL, CD, and crucially, the derivatives of CL and CD with respect to alpha, aspect ratio, and wing twist, back to the optimizer.

## 5. Expected Workflow (Conceptual)

The optimization process will generally follow these steps iteratively:

1.  **Initialization:** The optimizer starts with an initial guess for the design variables (alpha, aspect ratio, twist).
2.  **Aerodynamic Analysis:** OptVL takes these design variable values, updates the aircraft model and flight conditions, runs AVL, and computes the resulting CL and CD.
3.  **Constraint Evaluation:** The optimizer checks if the calculated CL meets the target (CL = 0.5).
    *   If CL is not 0.5, the optimizer notes this deviation.
4.  **Objective Evaluation:** The optimizer notes the current CD value.
5.  **Derivative-Informed Step:** Using the derivative information from OptVL, the optimizer determines how to change alpha, aspect ratio, and twist. The goal is to:
    *   Reduce CD.
    *   Drive CL closer to 0.5.
6.  **Iteration:** Steps 2-5 are repeated with the new set of design variables.
7.  **Convergence:** The process continues until the optimizer finds a set of design variables for which CD is minimized, and the CL constraint (CL = 0.5) is satisfied within a specified tolerance. The changes in design variables and objective function also become very small.

## 6. Python Script Structure (High-Level)

The Python script to implement this optimization will generally involve the following:

*   **Imports:** Import `OptVL`, the chosen optimizer library (e.g., `scipy.optimize`), and potentially `numpy` for numerical operations.
*   **Objective Function:** Define a Python function that:
    *   Accepts an array of design variables (alpha, AR, twist) as input.
    *   Uses OptVL to set these variables and run the aerodynamic analysis.
    *   Returns the drag coefficient (CD).
*   **Constraint Function:** Define a Python function that:
    *   Accepts an array of design variables as input.
    *   Uses OptVL to set these variables and run the aerodynamic analysis to get the lift coefficient (CL).
    *   Returns the difference between the actual CL and the target CL (e.g., `actual_CL - 0.5`). The optimizer will try to make this value zero.
*   **Optimizer Setup:**
    *   Provide the objective function and constraint function to the optimizer.
    *   Define bounds (min/max values) for each design variable.
    *   Set initial guess values for the design variables.
    *   Configure optimizer settings (e.g., algorithm, tolerances).
*   **Run Optimization:** Call the optimizer's main function to start the process.
*   **Results:** Print the optimal design variable values, the minimized CD, and the achieved CL. Optionally, plot the optimization history if available.

## 7. Placeholder for Code

The detailed Python code demonstrating this setup will be provided below.

### Python Implementation (using Scipy)
```python
[[PYTHON_CODE_PLACEHOLDER_SCIPY]]
```
*(Detailed explanation of the code will follow once it's inserted. This section will walk through how to set up the OptVL object, define the geometry and flight conditions, create the objective and constraint functions using OptVL's API, and then call the Scipy optimizer.)*

## 8. Interpreting Results

After the optimization run completes, you should examine:

*   **Optimal Design Variables:** The final values for angle of attack, aspect ratio, and wing twist. These represent the configuration that best meets the problem's objective and constraint.
*   **Minimized Drag Coefficient (CD):** The lowest CD value achieved.
*   **Achieved Lift Coefficient (CL):** This should be very close to the target CL of 0.5, confirming the constraint has been met.
*   **Optimization History (Optional):** If your optimizer provides it (many Scipy optimizers do via a callback function), plotting the history of how the design variables, CD, and CL changed over iterations can be very insightful. This helps visualize the convergence process.

## Conclusion

This tutorial outlines the conceptual framework for performing a lift-constrained drag minimization of a wing using OptVL and an external optimizer. It demonstrates how to define the problem, select design variables, and understand the roles of both OptVL and the optimizer in the workflow.

Once the Python code is implemented, you will be able to run this optimization yourself. We encourage you to then experiment by:
*   Changing the target CL.
*   Trying different design variables (if OptVL supports them easily, e.g., taper ratio).
*   Modifying the bounds of the existing design variables.
*   Using a different wing geometry as the baseline.

This hands-on approach will help solidify your understanding of how OptVL can be applied to practical aerodynamic design optimization problems.The file `docs/tutorial_lift_constrained_drag_min.md` has been successfully created with the detailed descriptive content for the tutorial.

Here's a summary of what was done:
1.  Confirmed the existence of `examples/rectangle.avl`, which will be used as the baseline geometry for the tutorial.
2.  Used the `create_file_with_block` tool to create `docs/tutorial_lift_constrained_drag_min.md`.
3.  Populated the file with the following sections as requested:
    *   Title: Tutorial: Lift-Constrained Drag Minimization of a Wing
    *   1. Problem Definition (Objective, Constraint, Aircraft using `examples/rectangle.avl`, Flight conditions)
    *   2. Design Variables (Angle of Attack, Wing Aspect Ratio, Linear Wing Twist with ranges/initial values)
    *   3. Optimization Setup Overview (Role of optimizer, OptVL providing aero analysis and derivatives)
    *   4. OptVL's Role (Loading geometry, setting conditions, applying variables, executing AVL, returning results)
    *   5. Expected Workflow (Conceptual iterative process)
    *   6. Python Script Structure (High-level points for imports, objective/constraint functions, optimizer setup, run, results)
    *   7. Placeholder for Code (`[[PYTHON_CODE_PLACEHOLDER_SCIPY]]` with a brief note)
    *   8. Interpreting Results (What to look for, optional plotting)
    *   Conclusion (Summary, encouragement for user modification)
4.  The content is focused on the conceptual setup, with the actual Python code to be added later.
5.  Markdown formatting has been used, and the language is aimed at the target audience.

The subtask is complete.
