# Optimizing with OpenMDAO

If you don't already have OpenMDAO, you can install it with `pip install OpenMDAO`

OpenMDAO doesn't actually do any of the optimization, it is simply a convenient framework for forming models for optimization. 
It installs SciPy along with its optimizer as a dependency during the installation process. 
If you want access to additional optimizers, I recommend the [pyOptSparse library](https://github.com/mdolab/pyoptsparse).


## Adding OptVL to a model 

OptVL comes with OpenMDAO components to make it simple to include it into an OpenMDAO model. 
All the hard/tedious work of specifying how the gradient functions should be called is already done for you in this group. 
Simply import the group and add it to your model. 
```python
import openmdao.api as om
from optvl import OVLGroup

model = om.Group()
model.add_subsystem("ovlsolver", OVLGroup(geom_file="aircraft.avl", mass_file="aircraft.mass"))
```
All the geometric variables and `alpha` and `beta` are automatically exposed as inputs to the group. 
The geometric variables of each surface are added as inputs of the form `<surface name>:<geo variable>`.
To turn on other inputs such as the reference values and case parameters, use the "input_param_vals" and "input_ref_vals" keywords to the group, respectively.
```python 
model.add_subsystem("ovlsolver", OVLGroup(geom_file="aircraft.avl",input_param_vals=True, input_ref_vals=True))
```
The force coefficients such as `CL`, `CD`, and `Cm` are automatically added as outputs.
To add stability derivatives or control surface derivatives as outputs, just use the `output_stability_derivs` and `output_con_surf_derivs` keyword arguments to the group to turn them on. 
```python
model.add_subsystem("ovlsolver", OVLGroup(geom_file="aircraft.avl", output_stability_derivs=True, output_con_surf_derivs=True))
```

## Output options 

To write out the AVL geometry and the TecPlot surface CP file at every time step, use the option `write_grid=True`. 
Additionally, if you are not interested in viewing the files in ParaView and instead only want to use TecPlot 360, then I also recommend adding `write_grid_sol_time=True`.

```python
model.add_subsystem("ovlsolver", OVLGroup(geom_file="aircraft.avl", write_grid=True))
# for easier postprocessing in TecPlot 360 (not compatible with ParaView)
model.add_subsystem("ovlsolver", OVLGroup(geom_file="aircraft.avl", write_grid=True, write_grid_sol_time=True ))

```

## Setting design variables
For simple optimization problems, the inputs and outputs can be directly added as design variables, objective functions, and constraints. 


```python
# directly setting the inputs as design variables
model.add_design_var("ovlsolver.Wing:aincs", lower=-10, upper=10)
model.add_design_var("ovlsolver.Elevator", lower=-10, upper=10)

# the outputs of OptVL can be used as constraints
model.add_constraint("ovlsolver.CL", equals=1.5)
model.add_constraint("ovlsolver.CM", equals=0.0, scaler=1e3)

# the scaler values bring the objective function to ~ order 1 for the optimizer
model.add_objective("ovlsolver.CD", scaler=1e3)
```

## Setting Optimizer options and running

After you have set up your model you will have to apply a driver. 
In our case our driver is the SLSQP optimizer of SciPy, but you could also apply other drivers to the model such as the DOE driver.

The next step is to tell OpenMDAO we are done with our problem, and it can be set up with `prob.setup(mode='rev')`

!!! Warning
    When using the OptVL group we must setup the problem with `mode='rev'`. 
    OptVL uses reverse mode derivatives because in general we will have more geometric variables than output functions of interest.

The line with `om.n2` creates an N2 diagram, which is helpful for examining our model. 
See this [page](https://openmdao.org/newdocs/versions/latest/features/model_visualization/n2_details/n2_details.html) for more information on the N2 diagram.


```python 

prob = om.Problem(model)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
prob.driver.options['tol'] = 1e-6
prob.driver.options['disp'] = True

prob.setup(mode='rev')
om.n2(prob, show_browser=False, outfile="vlm_opt.html")
prob.run_driver()

```


## Example script
Below is an OpenMDAO script analogous to the SciPy-based example. 


```python 
{% include "../examples/run_opt_om.py" %}
```
Which outputs

```
Optimization terminated successfully    (Exit mode 0)
            Current function value: 85.54669279688602
            Iterations: 19
            Function evaluations: 19
            Gradient evaluations: 19
Optimization Complete
-----------------------------------
ovlsolver.Elevator [0.10802368]
ovlsolver.Wing:aincs [11.46493221  1.92748591  2.14445549  2.27131259 -4.81944669]
ovlsolver.CD [0.08554669]
```
Which is about the same optimized result we got from the SciPy-only based optimization. 