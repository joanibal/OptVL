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



get_case_total_data -> get_total_forces()


get_case_coef_derivs -> get_control_stab_derivs 
get_case_stab_derivs -> get_stab_derivs

get_case_surface_data -> get_surface_forces
get_case_parameter -> get_parameter
set_case_parameter -> set_parameter
get_case_constraint -> get_constraint


get_strip_data -> get_strip_forces


add_constraint -> set_constraint
add_trim_condition -> set_trim_condition
add_constraint for variables

con options are
avl.f:      CONNAM(ICALFA) = 'alpha '
avl.f:      CONNAM(ICBETA) = 'beta  '
avl.f:      CONNAM(ICROTX) = 'pb/2V '
avl.f:      CONNAM(ICROTY) = 'qc/2V '
avl.f:      CONNAM(ICROTZ) = 'rb/2V '
avl.f:      CONNAM(ICCL  ) = 'CL    '
avl.f:      CONNAM(ICCY  ) = 'CY    '
avl.f:      CONNAM(ICMOMX) = 'Cl roll mom'
avl.f:      CONNAM(ICMOMY) = 'Cm pitchmom'
avl.f:      CONNAM(ICMOMZ) = 'Cn yaw  mom'
avl.f:        CONNAM(IC) = DNAME(N)
avl.f:         IF(INDEX(CONNAM(IC),CONN(1:NCONN)).NE.0) GO TO 25
avl.f:          WRITE(LU,1050) VARNAM(IV), CONNAM(IC), CONVAL(IC,IR)

- set_variable
ovl.add_variable("alpha", 0.00)

- set_constraint
ovl.set_constraint("Rudder", "Cn", 0.0)



removes executeRun. use execute_run instead

