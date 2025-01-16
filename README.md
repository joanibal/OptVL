# OptAVL
[![Downloads](https://static.pepy.tech/badge/pyavl-wrapper)](https://pepy.tech/project/pyavl-wrapper)
<!--- [![Downloads](https://static.pepy.tech/badge/pyavl-wrapper/month)](https://pepy.tech/project/pyavl-wrapper) --->
[Documentation](https://joanibal.github.io/OptVL/)

# AVL + python + optimization = OptVL

OptVL is a modified version of Mark Drela and Harold Youngren's famous AVL code with a python-wrapper and AD derivative routines for gradient-based optimization.
The python wrapper allows one to more easily conduct large parameter sweeps with AVL or to include AVL into a larger model. 
Together with the additionaly derivive routines, OptVL can be added to an optimzation loop for design refinement by a gradient-based optimizer. 
On top of this OptVL also has the ability to view the geometry and pressure distribution in python, Paraview, or Tecplot. 
<!-- Additionally, this wrapper provides access to more data than is available through traditional file output.  -->
<!-- Unlike in the output files which is limit to about 4 digits, the user has access to the full double precision data.  -->

## key features

- specify and modify the geometry from python 
- compute derivatives of total forces with respect to geometric variables
- derivatives of stability derivatives 
    - needed for neutral point constraint
    


# Installation
The best way to get OptVL is to install it through pip
```
pip install optvl
```
Windows, macOS, and Linux are all supported!


## Building locally
If you want to make OptVL locally then you have to clone the repository and use the following process.

In the root directory run
```
pip install . 
```

<!-- ## building step by step

To compile the avl library use 
```
make
```
This code has only been tested with gfortran and gnu95 compilers. 
If you want to use something besides gfortran you will have to modify the Makefile


and to install the pyavl package on your python path use 
```
pip install . 
```
or 
```
pip install . -e 
```
to install in development mode  -->

# Basic usage
The API of OptVL was made to mirror the usage of AVL through its text interface. 
The user loads in a geometry file, adds constraints, and then executes analysis runs.

The AVL wrapper is implemented in the `OVLSolver` class. 
To use this wrapper, first one must initialize the `AVLSolver` object with a geometry file and optionally a mass file. 
After, the user can add constraints and then execute the run to generate data. 
Below is a basic example of this workflow. 
```python
from optvl import AVLSolver
import numpy as np

ovl = OVLSolver(geo_file="aircraft.avl")
ovl.add_constraint("alpha", 0.00)

# control surface names from geometry file
ovl.add_constraint("Elevator", 0.00, con_var="Cm pitch moment")
ovl.add_constraint("Rudder", 0.00, con_var="Cn yaw moment")

ovl.set_case_parameter("Mach", 0.3)

# This is the method that acutally runs the analysis
ovl.execute_run()

print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    ovl.add_constraint("alpha", alpha)
    ovl.execute_run()
    run_data = ovl.get_case_total_data()
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )

print("----------------- CL sweep ----------------")
print("   Angle        Cl           Cd          Cdff          Cdv          Cm")
for cl in np.arange(0.6,1.6,0.1):
    ovl.add_trim_condition("CL", cl)
    ovl.execute_run()
    run_data = ovl.get_case_total_data()
    alpha = ovl.get_case_parameter("alpha")
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )
```

# Parameter sweep example

## taper ratio sweep


# Optimization example

## optimize twist distribution

