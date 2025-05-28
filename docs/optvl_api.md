# Comparison of Interface with AVL

The API was made to mirror that of AVL.
However, instead of using nested layers of options, users call the methods for getting/setting options and running analysis on `OVLSolver` objects.
The tables below relate AVL commands and the corresponding method in OptVL.
See the limitations sections for more information about commands that are not available in OptVL.

|action| AVL's  "AVL" level command| OptVL API call|
|-----|--|--|
|load a geometry |LOAD <geo file>|OVLSolver(geo_file=<geo file>)|
|load a mass file |MASS <mass file> |OVLSolver(geo_file=<geo file>, mass_file=<mass file>) |
|load a run case file| CASE <case file>| not supported|
| apply mass file | MSET 0 | done automatically |


<!-- 
The commands from the oper and mode menus are available 
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
``` -->


|action| AVL's "OPER" command| OptVL API call|
|-----|--|--|
|setting the angle of attack|a a <angle>| ovl.set_constraint("alpha", <angle>)|
| set variable such that constraint = val | <variable> <constraint> <val> | ovl.set_constraint(<variable>, <val>, con_var=<constraint>) |
| set CL  constraint|  c1; c 1.3| ovl.set_trim_condition("CL", 1.3)|
| run an analysis | x | ovl.execute_run() |
| after an analysis | FT |  ovl.get_total_forces() |
| get strip force data | ST | ovl.get_strip_forces() |
| get shear moment distribution | VM | ovl.get_strip_forces() |
| get control surface derivatives (e.g. dCL/dElevator)| ST | ovl.get_control_derivs() |
| get stability derivatives | ST | ovl.get_stab_derivs()|
| get stability derivatives in the body axis| SB | - |
| get/set reference data | RE | ovl.get_reference_data/set_reference_data()|
| get/set  design variables specified in AVL file | DE | -|
| get surface forces | FN | ovl.get_surface_forces() |
| get force distribution on strips| FS| ovl.get_strip_forces() |
| get force distribution on elements | FE | - |
| get body forces| FB | -|
| get high moments| HM | ovl.get_hinge_moments() |


|action| AVL's mode level command| OptVL API call|
|-----|--|--|
| get/set case parameters |M <var> <value>| ovl.get_reference_data(<var>)/set_reference_data(<var>, <value>)|
| get system matrix | S | ovl.get_system_matrix()|
| get eigenvalues| W | ovl.get_eigen_values()|
| get eigenvectors| not supported | ovl.get_eigen_vectors()|
| execute eigenmode calculation | N | ovl.execute_eigen_mode_calc() |

## Limitations
1. OptVL does not support multiple run cases since this would make the wrapping and derivative code more complex. Instead, create multiple solver instances and apply different parameters to each to replicate this functionality. 
2. There is no single precision version of OptVL available for download. You could compile one yourself if you really need this.
3. No support for working with design variables set in the AVL geometry file since we use a different system for modifying the geometry that allows the user to change any geometric parameter.
