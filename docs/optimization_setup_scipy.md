# Optimizing with SciPy

Before we can use SciPy we have to install it. 
Thankfully, SciPy is one of the best supported python packages in the world so to install it once can just use `pip install scipy`.

SciPy has many optimizers available, but I'm going to focus on using SLSQP since it can use the gradient information from OptVL and supports constraints. 
To use SciPy's SLSQP we will need to supply it with custom objective and constraint functions as well as the derivatives for each. 
These functions need to take in the design variables and apply them to our OptVL solver. 
The snippet below provides an example of an objective function.
```python
 def objective_function(x):
    ovl.set_constraint("Elevator", x[0])
    ovl.set_surface_params({"Wing":{"aincs":x[1:]}})
    
    ovl.execute_run()
    cd = ovl.get_total_forces()['CD']
    print(x, cd)

    return cd
```
Note, the objective function is specified by the return value of the function. 
If you wanted to save data about each iteration or write output the objective function would be a good place add that functionality. 
To supply the gradient information to SLSQP we have to define another function that returns the gradients for a given design variable vector.
```python
def objective_gradient(x):
    # Partial derivatives of the objective_function
    ovl.set_constraint("Elevator", x[0])
    ovl.set_surface_params({"Wing":{"aincs":x[1:]}})
    
    ovl.execute_run()
    
    sens = ovl.execute_run_sensitivies(['CD'])
    dcd_dele = sens['CD']['Elevator']
    dcd_daincs = sens['CD']['Wing']['aincs']
    
    # concatinate the two and return the derivs
    return np.concatenate(([dcd_dele], dcd_daincs))
```
The function `ovl.execute_run_sensitivies(['CD'])` does all the necessary work to compute the derivatives for the given list of functions. 
We just need to parse the `sens` dictionary for the derivatives with respect to the design variables we are interested in.

One also needs to repeat the process to create functions for the value and gradients of the constraints.

## Example script
The script below shows the full optimization script for minimizing the drag of the example aircraft in trim with respect to the wing twist and elevator position. 

```python 
{% include "../examples/run_opt_scipy.py" %}
```

Which should produce the output after the optimization
```
Optimization terminated successfully    (Exit mode 0)
            Current function value: 0.08551861616299436
            Iterations: 56
            Function evaluations: 58
            Gradient evaluations: 56
Optimization result:
Elevator: 0.09925034689068805
Wing twist: [11.51424344  2.00085521  2.20301866  2.32160102 -4.79748301]
Final CD: 0.08551861616299436
```