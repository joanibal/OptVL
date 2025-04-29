# Changes from pyAVL to OptVL

When I first started working on this project, I called it pyAVL. 
I changed the name to OptVL for two reasons
1. There were many other pyAVL's and people who were looking for it had trouble finding it
2. I added lots of optimization features and wanted the name to reflect this development
While I changed the name, I decided to also make some changes to the API to better match AVL's interface. 
The table below lists the changes. 


| pyAVL Method              | OptVL Equivalent Method   |
|-------------------------------|-------------------------------|
| `get_case_total_data`         | `get_total_forces()`          |
| `get_case_coef_derivs`        | `get_control_stab_derivs`     |
| `get_case_stab_derivs`        | `get_stab_derivs`             |
| `get_case_surface_data`       | `get_surface_forces`          |
| `get_case_parameter`          | `get_parameter`               |
| `set_case_parameter`          | `set_parameter`               |
| `get_case_constraint`         | `get_constraint`              |
| `get_strip_data`              | `get_strip_forces`            |
| `add_constraint`              | `set_constraint`              |
| `add_trim_condition`          | `set_trim_condition`          |
| `executeRun`                  | use `execute_run` instead     |

The output dictionaries for stability and control surface derivatives are also now flat. 
For example instead of `cs_derivs['CL']['Elevator']` the key is now `cs_derivs['dCL/dElevator']`
