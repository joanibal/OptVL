# Tutorial: Maximizing L/D Ratio with Geometric Variables (Scipy)

## 1. Introduction

The Lift-to-Drag (L/D) ratio is a critical metric in aircraft performance, directly influencing factors like range and endurance. A higher L/D ratio generally means better aerodynamic efficiency.

This tutorial demonstrates how to maximize the L/D ratio by simultaneously optimizing two key design variables: the angle of attack (`alpha`) and the wing aspect ratio (`AR`). We will use OptVL to perform the aerodynamic analyses and Scipy's Sequential Least SQuares Programming (SLSQP) optimizer to find the optimal design.

A core aspect of this tutorial is showcasing how to handle optimizations involving geometric modifications, specifically changes to aspect ratio, and how to compute the necessary sensitivities (gradients) for the optimizer, even when the conceptual design variable (AR) is not a direct input to the underlying analysis tool's sensitivity calculations.

## 2. Problem Setup

*   **Objective:** Maximize the lift-to-drag ratio (CL/CD).
*   **Design Variables:**
    *   **Angle of Attack (`alpha`):** Measured in degrees.
        *   Example Bounds: `(-5.0, 10.0)` degrees.
    *   **Wing Aspect Ratio (`AR`):** A dimensionless parameter (AR = span² / Sref).
        *   Example Bounds: `(4.0, 12.0)`.
*   **Aircraft:** We will use a standard AVL geometry file. The example script `examples/run_max_ld_ratio_scipy.py` is configured to use `aircraft.avl`. This file should be accessible to the script, typically located in an `examples/` directory.
*   **Flight conditions:** A representative Mach number is specified for the analysis. The example script uses Mach 0.1.
*   **Important Note on Aspect Ratio Modification:**
    Changing the aspect ratio of a wing while keeping its reference area (Sref) constant requires adjustments to both its span and chord dimensions. The example script `examples/run_max_ld_ratio_scipy.py` implements this by modifying specific sectional `SPAN` and `CHORD` parameters within OptVL. This approach is highly dependent on:
    1.  The specific structure of the input `aircraft.avl` file (e.g., how its sections are defined).
    2.  OptVL's Application Programming Interface (API) for geometry manipulation (i.e., how `set_parameter` calls translate to changes in the AVL model).
    The script includes detailed comments on these assumptions.

## 3. OptVL's Role

In this optimization task, OptVL is responsible for:

1.  **Loading Geometry:** Reading the baseline aircraft geometry from the `aircraft.avl` file.
2.  **Setting Flight Conditions:** Establishing the specified Mach number for the analysis.
3.  **Applying Design Variables:**
    *   Setting the current angle of attack (`alpha`).
    *   Modifying the wing geometry to achieve the target Aspect Ratio (`AR`). This involves:
        *   Calculating the new required span and mean geometric chord to maintain the aircraft's initial reference area (Sref).
        *   Using OptVL's `set_parameter` function to adjust underlying geometric properties (e.g., `SPAN` and `CHORD` of specific wing sections) to realize the new AR.
4.  **Executing Analysis:** Running the core AVL routines to perform the aerodynamic calculations for the current aircraft configuration and flight state, yielding total forces like CL and CD.
5.  **Providing Sensitivities:** Calculating and providing the sensitivities (derivatives) of CL and CD. This includes:
    *   Sensitivity with respect to `alpha`.
    *   Sensitivities with respect to the fundamental geometric parameters that were modified to achieve the AR change (e.g., sectional `SPAN` and `CHORD` values). These are then used to compute the sensitivity with respect to AR via the chain rule.

## 4. Optimization Workflow with Scipy

The optimization is driven by Scipy's `minimize` function, specifically using the SLSQP method, which is suitable for constrained optimization problems (though this example primarily uses bounds, SLSQP is robust).

*   **Objective Function:**
    A Python function is defined that:
    1.  Accepts an array `x` containing the current values of the design variables: `[alpha, AR]`.
    2.  Instructs OptVL to set these design variables (involving geometry modification for AR).
    3.  Commands OptVL to run the aerodynamic analysis and retrieve CL and CD.
    4.  Returns the value `-CL/CD`. Scipy's optimizers work by minimizing a function, so maximizing `CL/CD` is equivalent to minimizing `-CL/CD`.

*   **Gradient Function:**
    A separate Python function calculates the gradient of the objective function (`-CL/CD`) with respect to each design variable `[alpha, AR]`. This is crucial for the efficiency of gradient-based optimizers like SLSQP.
    *   **Gradient w.r.t. `alpha` (d(obj)/dalpha):** This derivative is typically obtained directly from OptVL's sensitivity analysis results for CL and CD with respect to alpha, and then applying the quotient rule for `-CL/CD`.
    *   **Gradient w.r.t. `AR` (d(obj)/dAR):** This is more complex because AR is not usually a direct input variable for which AVL (and thus OptVL) calculates sensitivities. Instead, AR is a conceptual variable that is realized by changing other fundamental geometric parameters (like sectional span and chord). The example script calculates `d(obj)/dAR` using the **chain rule**:
        1.  OptVL provides `dCL/d(geom_param)` and `dCD/d(geom_param)` where `geom_param` might be a sectional span or chord.
        2.  Mathematical expressions are derived for `d(geom_param)/dAR` (i.e., how much a specific sectional span or chord needs to change for a unit change in AR, while keeping Sref constant).
        3.  The chain rule combines these: `dCL/dAR = sum_over_geom_params( (dCL/d(geom_param)) * (d(geom_param)/dAR) )`. A similar calculation is done for `dCD/dAR`.
        4.  Finally, these are used with the quotient rule to find `d(-CL/CD)/dAR`.
        *This is an advanced but essential technique for optimizing designs using higher-level conceptual variables. The example script `examples/run_max_ld_ratio_scipy.py` contains detailed comments explaining its specific chain rule implementation.*

*   **Scipy `minimize` Call:**
    The script calls `scipy.optimize.minimize` configured with:
    *   The objective function.
    *   The gradient function (`jac=objective_gradient`).
    *   Bounds for each design variable (`alpha` and `AR`).
    *   An initial guess for the design variables.
    *   The chosen method (`'SLSQP'`).

## 5. Example Script

The Python script `examples/run_max_ld_ratio_scipy.py` implements this L/D maximization. Here is the content of the script:

```python
{% include "../examples/run_max_ld_ratio_scipy.py" %}
```

This script performs the following key operations:
*   Initializes the OptVL solver (`OVLSolver` or a placeholder if the library is not found).
*   Performs a baseline run to fetch initial geometric properties like Sref, span (b), and reference chord (Cref), which are crucial for consistent AR modification.
*   Defines the `objective_function` which applies design variables, runs OptVL, and returns `-L/D`.
*   Defines the `objective_gradient` function. This function includes the logic for:
    *   Directly using `alpha` sensitivities from OptVL.
    *   Implementing the chain rule to calculate sensitivities with respect to `AR` based on the sensitivities of underlying geometric parameters (sectional span and chord) that OptVL modifies.
*   Sets up initial guesses, bounds for `alpha` and `AR`.
*   Calls `scipy.optimize.minimize` using the `SLSQP` algorithm, providing both the objective and gradient functions.
*   Prints the optimization results.

## 6. Running the Tutorial & Interpreting Results

To run this tutorial, execute the script from your terminal (assuming you are in the main repository directory or have adjusted paths accordingly):

```bash
python examples/run_max_ld_ratio_scipy.py
```

**What to look for in the output:**

*   **Optimizer's Progress:** Scipy's `minimize` function (with `options={'disp': True}`) will print messages at each iteration, showing the current objective function value and other diagnostic information.
*   **Script Print Statements:** The example script includes print statements within the `objective_function` and `objective_gradient` functions. These will show:
    *   The current design variables (`[alpha, AR]`) being evaluated.
    *   The resulting CL, CD, and L/D ratio for those variables.
    *   The calculated gradient components `[d(obj)/dalpha, d(obj)/dAR]`.
    This detailed logging is very useful for monitoring the optimizer's behavior and debugging.
*   **Final Results:** Upon completion, the script will print:
    *   Whether the optimization was successful.
    *   The optimal values found for `alpha` and `AR`.
    *   The maximum L/D ratio achieved.
    *   A verification run showing CL, CD, and L/D using the optimal variables.

## 7. Key Learning Points & Considerations

*   **Importance of L/D:** Maximizing L/D is a fundamental task in aerodynamic design for improving aircraft efficiency.
*   **Geometric Optimization:** Optimizing geometric variables like Aspect Ratio often requires a deeper understanding of how these conceptual variables map to concrete parameters in the analysis tool (OptVL/AVL).
*   **Chain Rule for Gradients:** When an optimizer requires gradients with respect to a design variable (like AR) that is not a direct input for sensitivity analysis in the underlying tool, the chain rule becomes indispensable. This involves:
    1.  Identifying the fundamental parameters the tool *can* vary and provide sensitivities for (e.g., sectional chord, span).
    2.  Mathematically relating the conceptual variable to these fundamental parameters.
    3.  Applying the chain rule to combine these relationships with the sensitivities of the fundamental parameters.
    *The example script provides a concrete case of this. Refer to its comments for the detailed derivation and assumptions about how OptVL handles geometry and sensitivities.*
*   **Accuracy is Key:** The success of gradient-based optimization heavily relies on:
    *   The accuracy of the aerodynamic sensitivities provided by OptVL.
    *   The correctness of the chain rule implementation if used.
*   **OptVL Installation:** The example script includes a placeholder `OVLSolverPlaceholder`. For actual use, ensure that the OptVL library is correctly installed and that `OVLSolver` (or the relevant class name) is imported properly. The script attempts to import the real solver first.
*   **AVL File Structure:** The way geometric parameters (like sectional span and chord) are modified in the script is tied to assumptions about the structure of the `aircraft.avl` file. Using a different AVL file might require changes to the geometry modification logic in the script.

## Conclusion

This tutorial demonstrated a practical example of aircraft aerodynamic optimization: maximizing the L/D ratio by varying angle of attack and wing aspect ratio. It highlighted the role of OptVL in providing aerodynamic data and sensitivities, and how Scipy can be used to drive the optimization.

A key takeaway is the methodology for handling geometric design variables like aspect ratio, particularly the use of the chain rule for gradient calculations when direct sensitivities for such variables are not available.

We encourage users to:
*   Study the `examples/run_max_ld_ratio_scipy.py` script in detail, paying close attention to the `apply_design_variables` and `objective_gradient` functions.
*   Experiment with different bounds, initial guesses for `alpha` and `AR`.
*   Consider adapting the script for a different AVL file, but be mindful that this will likely require adjustments to the geometry modification logic (i.e., which `set_parameter` calls are made and how their values are calculated from AR).
Happy optimizing!
