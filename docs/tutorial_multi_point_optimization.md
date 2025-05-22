# Tutorial: Multi-Point Aerodynamic Optimization with OpenMDAO

## 1. Introduction

Aircraft are rarely designed to perform optimally at only a single flight condition. Instead, they often need to operate efficiently and effectively across a range of scenarios, such as cruise, climb, loiter, or high-speed maneuvers. Designing for this versatility is a complex challenge.

Multi-point optimization (MPO) is an engineering technique used to find a design that offers the best possible compromise across these varied operating conditions. Instead of optimizing for a single point, MPO considers performance at multiple points simultaneously, often by minimizing a weighted sum of objectives or by satisfying constraints at each point.

This tutorial demonstrates how to perform a multi-point aerodynamic optimization using OptVL and OpenMDAO. Specifically, we will optimize a shared set of wing twist angles (`wing_aincs`) and individual elevator settings for two distinct flight conditions (cruise and loiter). The goal is to minimize a weighted average of the drag coefficients from these two conditions, while ensuring that a target lift coefficient (CL) is achieved at each condition.

## 2. Problem Setup

*   **Objective:** Minimize a weighted sum of drag coefficients from two flight conditions:
    `weighted_CD = w1 * CD_cruise + w2 * CD_loiter`
    where `w1` and `w2` are weighting factors representing the relative importance or time spent in each condition.

*   **Design Variables:**
    *   **Shared Design Variable:**
        *   `wing_aincs`: A set of wing incidence angles (representing wing twist) distributed along the span. These twist angles are part of the aircraft's physical design and are therefore common to both flight conditions.
    *   **Point-Specific Design Variables:**
        *   `Elevator_cruise`: Elevator deflection angle (in degrees) for the cruise flight condition. This is used to trim the aircraft or achieve the target CL for cruise.
        *   `Elevator_loiter`: Elevator deflection angle (in degrees) for the loiter flight condition. This is used to trim the aircraft or achieve the target CL for loiter.

*   **Constraints:**
    *   Achieve the target lift coefficient for Condition 1 (Cruise): `cruise_point.CL = target_CL_cruise`.
    *   Achieve the target lift coefficient for Condition 2 (Loiter): `loiter_point.CL = target_CL_loiter`.

*   **Flight Conditions:**
    *   **Condition 1 (Cruise):**
        *   Mach Number: `Mach_cruise` (e.g., 0.8)
        *   Target Lift Coefficient: `target_CL_cruise` (e.g., 0.5)
        *   Weighting Factor: `w1` (e.g., 0.7)
    *   **Condition 2 (Loiter):**
        *   Mach Number: `Mach_loiter` (e.g., 0.4)
        *   Target Lift Coefficient: `target_CL_loiter` (e.g., 0.7)
        *   Weighting Factor: `w2` (e.g., 0.3)

*   **Aircraft Model:** The optimization will use a geometry file named `aircraft.avl`. OptVL will use this file for its aerodynamic analyses.

## 3. OpenMDAO Model Structure

OpenMDAO (Open Multi-Disciplinary Analysis and Optimization) is exceptionally well-suited for MPO problems due to its modular, component-based architecture. We can construct a model where each flight condition is analyzed by a dedicated component, and their outputs are then combined to form the overall objective.

The OpenMDAO model for this tutorial is structured as follows:

*   **Main Group (`MultiPointOptimizationGroup`):** This is the top-level group that encapsulates the entire multi-point problem. It contains:
    *   **Two instances of OptVL's `OVLGroup`:**
        *   `cruise_point`: An `OVLGroup` instance configured for Condition 1 (Cruise).
        *   `loiter_point`: An `OVLGroup` instance configured for Condition 2 (Loiter).
        Each `OVLGroup` is responsible for running an OptVL analysis for its specific flight condition.
    *   **`IndepVarComp` (Independent Variable Components):**
        *   A `shared_ivc` is used to define the source for shared design variables, specifically `wing_aincs_source`. This output is then connected to both `OVLGroup` instances.
        *   A `point_ivc` is used to define point-specific independent variables. These include the Mach numbers for each condition (`Mach_cruise`, `Mach_loiter`) and the initial guesses for the elevator settings (`Elevator_cruise`, `Elevator_loiter`). The Mach numbers are fixed parameters for this problem, while the elevator settings are design variables.
    *   **Custom `ExplicitComponent` (`WeightedDragObjective`):**
        *   This component is created to calculate the weighted average drag. It takes the `CD` outputs from `cruise_point` and `loiter_point` as inputs and computes `avg_CD = w1*CD_cruise + w2*CD_loiter`.

*   **Key Connections within the Model:**
    *   The shared design variable `shared_ivc.wing_aincs_source` is connected to the `Wing:aincs` input of both `cruise_point` and `loiter_point`. This ensures both analyses use the same wing twist.
    *   The Mach number outputs from `point_ivc` (`Mach_cruise`, `Mach_loiter`) are connected to the respective `Mach` inputs of `cruise_point` and `loiter_point`.
    *   The elevator setting outputs from `point_ivc` (`Elevator_cruise`, `Elevator_loiter`) are connected to the respective `Elevator` inputs of `cruise_point` and `loiter_point`. (Note: These are also designated as design variables for the optimizer).
    *   The drag coefficient outputs (`CD`) from `cruise_point` and `loiter_point` are connected to the `CD_cruise` and `CD_loiter` inputs of the `WeightedDragObjective` component.

*   **Defining Design Variables, Objective, and Constraints in OpenMDAO:**
    Once the model structure is built, the optimization problem is defined at the `Problem` level:
    *   `prob.model.add_design_var('shared_ivc.wing_aincs_source', ...)`: Tells OpenMDAO to treat the `wing_aincs_source` output of the `shared_ivc` as a design variable that the optimizer can change. Similar calls are made for `point_ivc.Elevator_cruise` and `point_ivc.Elevator_loiter`.
    *   `prob.model.add_objective('objective_comp.avg_CD')`: Sets the `avg_CD` output of the `WeightedDragObjective` component as the quantity to be minimized.
    *   `prob.model.add_constraint('cruise_point.CL', equals=target_CL_cruise)`: Constrains the `CL` output of the `cruise_point` analysis to be equal to the target CL for that condition. A similar constraint is added for `loiter_point.CL`.

## 4. OptVL's Role (`OVLGroup`)

Each instance of `OVLGroup` (provided by OptVL for OpenMDAO integration) acts as a black box for one specific aerodynamic analysis point:

*   It is configured with the aircraft geometry (`aircraft.avl`) and a unique name for the flight condition it represents.
*   It receives inputs such as:
    *   `Mach`: Mach number for its specific flight condition.
    *   `Elevator`: Elevator deflection angle (point-specific design variable).
    *   `Wing:aincs`: The array of wing incidence/twist angles (shared design variable).
*   Based on these inputs, it runs an OptVL/AVL analysis.
*   It then outputs the resulting aerodynamic coefficients like `CL`, `CD`, `CM`, and potentially other data for its configured flight condition.
*   Crucially for gradient-based optimization, `OVLGroup` is also responsible for computing and providing the partial derivatives (sensitivities) of its outputs (CL, CD, etc.) with respect to its inputs (Mach, Elevator, Wing:aincs, etc.). OpenMDAO uses these derivatives to efficiently find the optimal solution.

## 5. Example Script

The Python script `examples/run_opt_multipoint_om.py` implements this multi-point optimization. Here is the content of the script:

```python
{% include "../examples/run_opt_multipoint_om.py" %}
```

This script performs the following key operations:
*   Defines the flight conditions (Mach numbers, target CLs, objective weights).
*   Defines the `WeightedDragObjective` OpenMDAO component.
*   Constructs the `MultiPointOptimizationGroup` by:
    *   Instantiating two `OVLGroup` components (`cruise_point`, `loiter_point`).
    *   Setting up `IndepVarComp`s for shared and point-specific variables.
    *   Connecting all the inputs and outputs between these components as described above.
*   Configures the OpenMDAO `Problem` by:
    *   Adding the design variables (`wing_aincs_source`, `Elevator_cruise`, `Elevator_loiter`).
    *   Setting the objective (`objective_comp.avg_CD`).
    *   Adding the CL constraints for both flight points.
*   Initializes and runs the `ScipyOptimizeDriver` (using SLSQP).
*   Prints the results of the optimization.

## 6. Running the Tutorial & Interpreting Results

To run this tutorial, navigate to the directory containing the repository and execute the script from your terminal:

```bash
python examples/run_opt_multipoint_om.py
```

**What to look for in the output:**

*   **Optimizer's Progress:** If `disp=True` is set for the Scipy driver, OpenMDAO will print information at each iteration of the optimization, including function calls, objective value, and constraint satisfaction.
*   **Final Optimal Values:** Once the optimization converges, the script will print:
    *   The optimal values for the shared wing twist angles (`wing_aincs_source`).
    *   The optimal elevator deflection for the cruise condition (`Elevator_cruise`).
    *   The optimal elevator deflection for the loiter condition (`Elevator_loiter`).
*   **Minimized Objective:** The final value of the `weighted_CD` (the objective function).
*   **Constraint Satisfaction:** The script also prints the final CL values for both the cruise and loiter points, allowing you to verify that they meet their respective target CLs.

## 7. Key Learning Points

*   **MPO for Versatility:** Multi-point optimization is a powerful and often necessary approach for designing aircraft that must perform well across a spectrum of flight conditions.
*   **OpenMDAO for MPO:** OpenMDAO's component-based architecture and its ability to manage complex data flows make it an excellent framework for setting up and solving MPO problems. You can clearly separate the analysis for each flight point while managing shared design parameters.
*   **Shared vs. Point-Specific Variables:**
    *   Shared design variables (like the physical wing twist) influence performance in all considered flight conditions. The optimizer must find a single set of these variables that contributes best to the overall multi-point objective.
    *   Point-specific design variables (like elevator deflections or throttle settings) can be adjusted independently for each flight condition to meet specific targets (e.g., trimming the aircraft to a target CL).
*   **Reusability of `OVLGroup`:** The `OVLGroup` component from OptVL (or its equivalent) can be instantiated multiple times within the same OpenMDAO model, each instance configured for a different analysis point or flight condition, using the same underlying geometry file.
*   **Placeholder for OptVL:** The example script `run_opt_multipoint_om.py` is designed to be runnable even if OptVL is not fully installed, by using a `OVLGroupPlaceholder`. For actual aerodynamic optimization, a working OptVL installation and the correct import of `OVLGroup` are necessary.

## Conclusion

This tutorial has outlined the process of performing a multi-point aerodynamic optimization using OptVL and OpenMDAO. By considering multiple flight conditions simultaneously, engineers can arrive at designs that are more robust and versatile, leading to better overall aircraft performance.

OptVL, when integrated into OpenMDAO via components like `OVLGroup`, provides the aerodynamic analysis capabilities and crucial derivative information needed for efficient gradient-based MPO. OpenMDAO, in turn, offers the framework to orchestrate these multiple analyses and drive the design towards an optimal compromise.

We encourage users to explore this example further by:
*   Modifying the weighting factors (`WEIGHT_CRUISE`, `WEIGHT_LOITER`) to see how they influence the optimal design.
*   Changing the target CL values or Mach numbers.
*   Expanding the problem to include more flight conditions.
*   Introducing other shared or point-specific design variables that OptVL can control.

Happy optimizing!
