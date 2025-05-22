# Why OptVL? A Guide for AVL Users

## Introduction

If you're an aerospace engineer, you're likely familiar with Athena Vortex Lattice (AVL) and its strengths in aerodynamic analysis. AVL is a powerful tool, and OptVL is designed to extend its capabilities by bringing them into the Python ecosystem. The primary goal of OptVL is to leverage AVL's validated physics for broader applications, especially in areas like design optimization and extensive parametric studies. Think of OptVL as a bridge that connects AVL's core functionalities with the flexibility and power of Python.

## Key Advantages of OptVL

OptVL offers several advantages over using AVL in its traditional, standalone mode:

*   **Automation and Scripting:**
    *   OptVL allows you to write Python scripts to automate your AVL workflows. This means you can programmatically set parameters (like angle of attack, control surface deflections), run analysis cases, and extract results without manual intervention.
    *   This is a significant improvement over the manual input and output process of traditional AVL, which can be time-consuming and error-prone for large numbers of cases.
    *   This is particularly beneficial for tasks like generating aerodynamic databases, where hundreds or even thousands of AVL runs might be needed.

*   **Design Optimization:**
    *   OptVL is specifically built with optimization in mind. It's not just about running AVL cases, but about enabling design improvements.
    *   A key feature is its ability to calculate the derivatives (sensitivities) of aerodynamic forces and moments with respect to changes in design variables (e.g., wing span, sweep, airfoil parameters). These derivatives are crucial for efficient gradient-based optimization algorithms.
    *   OptVL is designed to integrate with popular Python optimization libraries such as Scipy and OpenMDAO, allowing you to use sophisticated optimizers to find the best designs. (More details on this integration can be found in other documentation pages).

*   **Integration with Python Ecosystem:**
    *   By bringing AVL into Python, OptVL allows you to connect AVL's analysis capabilities with the vast array of other Python libraries.
    *   For example, you can use `matplotlib` for advanced and customized plotting of aerodynamic data, `numpy` and `scipy` for complex numerical calculations beyond AVL's standard output, or even integrate with machine learning tools for surrogate modeling or data analysis.

*   **Enhanced Data Access:**
    *   OptVL can provide more direct and detailed access to internal AVL data structures and results than what is typically available through standard AVL output files (like `.st` or `.fs` files). This allows for more in-depth analysis and post-processing.

## Workflow Comparison (High-Level)

While OptVL introduces new capabilities, it aims to keep the core AVL concepts familiar. Here’s a high-level comparison of how common tasks are performed:

*   **Geometry Definition:**
    *   **AVL:** Geometry is defined in an AVL input file (e.g., `my_aircraft.avl`).
    *   **OptVL:** OptVL typically loads these existing AVL geometry files. The way you define your aircraft geometry remains largely the same.

*   **Defining an Analysis Case:**
    *   **AVL:** You manually enter commands in the AVL interface, such as `OPER` to enter the operating conditions menu, then set parameters like angle of attack (`A A 5.0`), sideslip angle, control surface deflections, etc.
    *   **OptVL:** You use Python functions to set these parameters. For example, `optvl_instance.set_parameter("alpha", 5.0)`.

*   **Running an Analysis:**
    *   **AVL:** You use the `X` (execute run) command.
    *   **OptVL:** You call a Python method, typically something like `results = optvl_instance.execute_run()`.

*   **Accessing Results:**
    *   **AVL:** You view results as text output in the console, or parse text files like `.st` (stability derivatives), `.fs` (force summary), or total forces files.
    *   **OptVL:** Results are returned as Python objects, such as dictionaries or NumPy arrays. This makes it much easier to access specific values, plot data, or perform further calculations programmatically. For instance, `CL = results['forces']['CL']`.

*   **Conceptual Mapping:**
    *   The core concepts and commands in AVL (like setting constraints, defining control surfaces, running different analysis types) have corresponding Python functions or methods in OptVL. This makes the transition more intuitive for experienced AVL users. For a detailed mapping of AVL commands to OptVL methods, please refer to the `optvl_api.md` document.

## How OptVL Extends AVL

OptVL isn't just a wrapper; it extends AVL's capabilities:

*   **Derivatives:** As mentioned, OptVL incorporates the calculation of aerodynamic derivatives using advanced methods like the adjoint method. This is a significant addition that makes gradient-based optimization feasible and efficient.
*   **Python Interface:** OptVL uses `f2py` (Fortran to Python interface generator) to wrap and call modified AVL Fortran routines directly. This allows for a tight integration and efficient communication between Python and the core AVL code.

## Conclusion

OptVL aims to empower you, the AVL user, by making AVL's robust aerodynamic analysis capabilities more flexible, automatable, and powerful through integration with the Python ecosystem. It opens the door to advanced design optimization, large-scale parametric studies, and seamless integration with other engineering tools.

We encourage you to explore the installation guide and basic usage tutorials next to get started with OptVL. For those interested in the direct mapping of AVL commands to OptVL functions, the `optvl_api.md` document provides a comprehensive reference.
Welcome to the enhanced capabilities of AVL with OptVL!
