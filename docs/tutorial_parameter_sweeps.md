# Tutorial: Design Space Exploration with Parameter Sweeps using OptVL

## 1. Introduction

Before embarking on formal optimization, it's often crucial to gain a deeper understanding of how different design parameters affect an aircraft's performance. This process is often referred to as "design space exploration." Parameter sweeps are a common and effective method for this. They involve systematically varying one or more design parameters across a defined range and observing their impact on key performance metrics.

This tutorial will guide you through performing a 2D parameter sweep. Specifically, we will vary the Angle of Attack (`alpha`) and the wing Aspect Ratio (`AR`) for an aircraft model. We'll use OptVL to automate the aerodynamic analyses for each combination of these parameters and Matplotlib to visualize the results, helping us understand the aerodynamic characteristics of the aircraft.

## 2. Purpose of Parameter Sweeps

Parameter sweeps serve several important purposes in the engineering design process:

*   **Gaining Design Sensitivities:** Understand how much key performance metrics (like CL, CD, L/D) change in response to changes in design parameters.
*   **Identifying Interesting Regions:** Pinpoint areas in the design space that exhibit desirable characteristics (e.g., high L/D, linear lift behavior) or problematic behavior (e.g., stall, high drag).
*   **Visualizing Trade-offs:** See how improving one aspect of performance might affect another when parameters are changed. For example, how does increasing aspect ratio affect both CL_max and CD_min?
*   **Generating Data for Surrogate Models:** Although this tutorial focuses on direct visualization, the data generated from extensive sweeps can be used to train surrogate models (simplified mathematical approximations of the complex analysis code), which can speed up later optimization tasks.
*   **Sanity Checking Models and Assumptions:** Unexpected results from a parameter sweep can often highlight issues with the analysis model, its setup, or underlying assumptions.

## 3. Setup for the Sweep

To conduct our parameter sweep, we need to define the parameters we'll vary and how we'll analyze their effects.

*   **Parameters to Sweep:**
    *   **Angle of Attack (`alpha`):** This will be swept across a range of values, for example, from -2 degrees to 10 degrees, to capture different flight attitudes.
    *   **Wing Aspect Ratio (`AR`):** We will use a set of discrete values for aspect ratio, such as 6, 8, 10, and 12. Aspect ratio (AR = span² / Sref) is a key geometric parameter influencing aerodynamic efficiency.

*   **Aerodynamic Analysis:** For each combination of `alpha` and `AR`, OptVL will be used to perform an aerodynamic analysis, calculating the lift coefficient (CL) and drag coefficient (CD). From these, the lift-to-drag ratio (L/D = CL/CD) will be derived.

*   **Aircraft Model:** The example script `examples/run_parameter_sweep.py` is configured to use the `aircraft.avl` geometry file. This file must be accessible to the script.

*   **Geometry Modification for AR:**
    Changing the aspect ratio of a wing typically involves modifying its span and/or chord dimensions. In this tutorial (and the accompanying script), we aim to change AR while keeping the wing's reference area (Sref) constant. This requires calculating the new span and mean geometric chord for each target AR value and then instructing OptVL to apply these changes to the AVL model. The exact OptVL calls (`set_parameter` for specific sectional `SPAN` and `CHORD` values) are dependent on the structure of the `aircraft.avl` file and OptVL's API for geometry manipulation. The example script includes detailed comments on these assumptions.

## 4. Script Workflow

The Python script (`examples/run_parameter_sweep.py`) automates the parameter sweep process with the following general workflow:

1.  **Initialization:**
    *   Import necessary libraries: `optvl.solver` (or a placeholder), `numpy` for numerical operations, and `matplotlib.pyplot` for plotting.
    *   Initialize an instance of the OptVL solver (`OVLSolver`).

2.  **Initial Geometry Properties:**
    *   Before starting the sweep, the script loads the `aircraft.avl` file and runs a baseline analysis.
    *   From this baseline, it retrieves and stores initial geometric properties like the reference area (Sref), initial span (b_initial), and initial reference chord (Cref_initial). These initial values are crucial for consistently applying changes to aspect ratio while maintaining Sref.

3.  **Parameter Ranges:**
    *   Define the range of values for `alpha` (e.g., using `numpy.linspace`).
    *   Define the list of discrete values for `AR`.

4.  **Nested Loops for Sweeping:**
    *   The script uses nested loops to iterate through every combination of the defined `AR` and `alpha` values.
    *   The outer loop iterates through each specified `AR` value.
    *   The inner loop iterates through each `alpha` value in its defined range.

5.  **Inside the Loops (for each parameter combination):**
    *   **Apply Aspect Ratio:** For the current `AR` from the outer loop, the script calculates the necessary new span and chord (keeping Sref constant) and uses OptVL's `set_parameter` calls to modify the geometry of the aircraft model.
    *   **Set Angle of Attack:** The current `alpha` value from the inner loop is set using OptVL's `set_parameter` method.
    *   **Run Analysis:** OptVL's `execute_run` method is called to perform the aerodynamic analysis.
    *   **Store Results:** The calculated CL, CD, and derived L/D ratio for the current `alpha` and `AR` combination are stored, typically in a list of dictionaries or a similar data structure. Progress is printed to the console.

6.  **Data Aggregation:** All results (alpha, AR, CL, CD, L/D for each run) are collected in a comprehensive data structure.

7.  **Plotting:** After the nested loops complete and all data points are gathered, Matplotlib is used to generate various plots to visualize the relationships between the parameters and the performance metrics.

## 5. Example Script

The Python script `examples/run_parameter_sweep.py` implements this parameter sweep. Here is the content of the script:

```python
{% include "../examples/run_parameter_sweep.py" %}
```

This script performs the key operations described in the workflow:
*   It initializes OptVL (or a placeholder for basic execution).
*   It defines the ranges for `alpha` and `AR`.
*   It uses nested loops to cover all combinations of these parameters.
*   Within the loops, it modifies the aircraft geometry to achieve the target `AR` (while keeping Sref constant) and sets the `alpha`.
*   It calls OptVL to run the analysis and collects CL, CD, and L/D.
*   Finally, it uses Matplotlib to generate and display several informative plots.

## 6. Visualizing Results (Types of Plots Generated by the Script)

The example script `examples/run_parameter_sweep.py` generates the following plots to help visualize the results of the parameter sweep:

*   **CL vs. Alpha (for each AR):**
    *   This plot shows the lift curve (CL as a function of `alpha`) for each aspect ratio tested.
    *   It helps visualize how the lift coefficient changes with angle of attack and how the lift curve slope and CL_max might be affected by aspect ratio.

*   **L/D vs. Alpha (for each AR):**
    *   This plot illustrates the lift-to-drag ratio (L/D) as a function of `alpha` for each AR.
    *   It is crucial for identifying the optimal angle of attack that yields the maximum L/D for each wing geometry. It also shows how the peak L/D value and the range of efficient operation change with AR.

*   **Drag Polars (CL vs. CD for each AR):**
    *   These are classic aerodynamic plots showing lift coefficient (CL) versus drag coefficient (CD). Each line represents a different aspect ratio.
    *   Drag polars provide a concise view of the overall aerodynamic efficiency of the aircraft across its lift range for different geometries. The "bottom" of each curve typically indicates the minimum drag for that configuration.

*   **Max L/D vs. AR:**
    *   This summary plot shows the maximum L/D ratio achieved for each aspect ratio tested.
    *   It clearly highlights how aspect ratio influences the peak aerodynamic efficiency of the wing, providing direct insight into the benefits or drawbacks of changing AR.

## 7. Key Learning Points

*   **Fundamental Exploration Tool:** Parameter sweeps are a foundational technique for understanding complex engineering systems and design trade-offs before (or alongside) formal optimization.
*   **Automation with OptVL:** OptVL can be easily scripted to automate the process of running numerous aerodynamic analyses for many different parameter combinations, making extensive sweeps feasible.
*   **Visualization for Insight:** Libraries like Matplotlib are essential for transforming raw data from sweeps into understandable visual plots, which are crucial for extracting insights and making informed design decisions.
*   **Careful Geometry Modification:** When sweeping geometric parameters like Aspect Ratio, it's vital to ensure that the modifications are applied correctly and consistently. The script's approach (modifying sectional span and chord to maintain Sref) and its underlying assumptions about the AVL file structure and OptVL's API should be well understood.
*   **Placeholder for Accessibility:** The example script includes a placeholder for `OVLSolver`. This allows users to run the script and see the data flow and plotting logic even without a full OptVL installation. However, for generating real aerodynamic data, a functioning OptVL environment is required.

## Conclusion

Parameter sweeps, when combined with tools like OptVL for analysis and Matplotlib for visualization, provide a powerful means to explore the design space and gain valuable insights into aircraft aerodynamic characteristics. This understanding is critical for making informed design choices and for setting up more targeted optimization studies.

We encourage users to:
*   Experiment with the `examples/run_parameter_sweep.py` script.
*   Modify the ranges for `alpha` and `AR`.
*   Adapt the script to sweep other geometric or flight condition parameters that OptVL can control.
*   Explore different visualization techniques or plot additional derived metrics relevant to their specific interests.

By systematically exploring the design space, engineers can build a strong intuition for their designs and lay a solid foundation for subsequent optimization efforts.
