# Basic Usage Overview

`OptVL` offers a simple and intuitive API, mirroring AVL's text interface.
The typical workflow involves loading a geometry file, adding constraints, and executing analysis runs.


The commands from the oper and mode menus are available 


```
  C1  set level or banked  horizontal flight constraints
  C2  set steady pitch rate (looping) flight constraints
  M odify parameters                                    

 "#" select  run case          L ist defined run cases   
  +  add new run case          S ave run cases to file   
  -  delete  run case          F etch run cases from file
  N ame current run case       W rite forces to file     

 eX ecute run case             I nitialize variables     

  G eometry plot               T refftz Plane plot       

  ST  stability derivatives    FT  total   forces        
  SB  body-axis derivatives    FN  surface forces        
  RE  reference quantities     FS  strip   forces        
  DE  design changes           FE  element forces        
  O ptions                     FB  body forces           
                               HM  hinge moments         
                               VM  strip shear,moment    
  MRF  machine-readable format CPOM OML surface pressures
```

|action| avl oper command| OptVL api call|
|-----|--|--|
|setting the angle of attack|a a 5|ovl.set_constraint("alpha", 5.0)|
| add that variable is set such that constraint = val | <variable> <constraint> <val> | 
| set CL  constraint|  c1; c 1.3| ovl.set_trim_condition("CL", 1.3)|
| run an analysis | x | ovl.execute_run() |
| after an analysis | FT |  get_total_forces() |
| get strip force data | ST | get_strip_forces() |
| get shear mommen distribution | VM | get_strip_forces() |
| get control surfaces derivatives (e.g. dCL/dElevator)| ST | get_control_derivs |
| get stability derivatives | ST | get_stab_derivs|
| get stability derivatives in the body axis| SB | - |
| get/set reference data | RE | get_reference_data/set_reference_data|
| get/set  design variables specified in AVL file | DE | -|
| look at options?? | O | - |
| get surface forces | FN | get_surface_forces |
| get force distribution on strips| FS| get_strip_forces |
| get force distribution on elements | FE | - |
| get body forces| FB | -|
| get high moments| HM | get_hinge_moments |



## Initializing and Setting up AVL Solver
To begin with `OptVL`, start by initializing the `AVLSolver` class:

```python
from optvl import OVLSolver
```
ovl = OVLSolver(geo_file="aircraft.avl")

After initializing, you can set up various constraints and parameters:

```python
ovl.set_constraint("alpha", 0.00)
ovl.set_constraint("Elevator", 0.00, con_var="Cm pitch moment")
ovl.set_constraint("Rudder", 0.00, con_var="Cn yaw moment")
ovl.set_parameter("Mach", 0.3)
```
## Running Analysis

Once you've set up the solver, running the analysis is straightforward:

```python
ovl.execute_run()
```

For a more detailed example and advanced use cases, see the analysis guide.

## Vizulization
