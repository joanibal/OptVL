# Analysis with OptVL

After initializing and setting up your `AVLSolver`, you can perform analysis tasks such as alpha sweeps and CL sweeps.

## Alpha Sweep

Here's an example of an alpha sweep:

```python
print("----------------- alpha sweep ----------------")
print("   Angle        Cl           Cd          Cdi          Cdv          Cm")
for alpha in range(10):
    avl_solver.set_constraint("alpha", alpha)
    avl_solver.execute_run()
    run_data = avl_solver.get_total_forces()
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )
```
## CL Sweep

Here's how you can perform a CL sweep:

```python
print("----------------- CL sweep ----------------")
print("   Angle        Cl           Cd          Cdff          Cdv          Cm")
for cl in np.arange(0.6,1.6,0.1):
    avl_solver.set_trim_condition("CL", cl)
    avl_solver.execute_run()
    run_data = avl_solver.get_total_forces()
    alpha = avl_solver.get_parameter("alpha")
    print(
        f' {alpha:10.6f}   {run_data["CL"]:10.6f}   {run_data["CD"]:10.6f}   {run_data["CDi"]:10.6f}   {run_data["CDv"]:10.6f}   {run_data["CM"]:10.6f}'
    )
```
