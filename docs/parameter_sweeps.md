# Modifying Geometry for Parameter Sweeps in OptVL

OptVL not only offers aerodynamic analysis capabilities but also provides tools to access and modify geometric parameters of your aircraft models. This documentation illustrates how users can retrieve and set geometry parameters for parameter sweeps.


## Setting geometry parameters

To apply geometry changes to the AVLSolver object, use the `set_surface_params` method.
```python
avl_solver.set_surface_params(data)
```
This method sets the surface parameters of the AVL model based on the provided dictionary, data.
The data pasted to this method must be a dictionary of surface names whos values are a dictionary of surface and section keywords.
An example of such a dictionary is
```python
data = {
    'Wing': {
        'angle': 3.5,  # a surface keyword example
        "chords": np.array([0.5, 0.4, 0.3, 0.2, 0.1]) #section keyword example
    }
}
```

## Getting geometry parameters

To retrieve surface parameters from your AVL model, use the `get_surface_params` method.
By default, this method only returns data about the geometry of the surface, but information about the paneling and control surfaces can also be included by passing the corresponding flags like so. 
```python
surf_data = avl_solver.get_surface_params(
    include_geom=True,
    include_panneling=True,
    include_con_surf=True
)
```
The data from `get_surface_params` come directly from the geometry file used by AVL.
See the [AVL user guide](https://web.mit.edu/drela/Public/web/avl/avl_doc.txt) for more information about all the possible variables.
For most use cases, it is probable that you will only need to interact with the geometric variables below. 

| Variable    | Description                       |
|-------------|-----------------------------------|
| scale       | Scale factor.                     |
| translate   | Translation vector.               |
| angle       | Surface angle.                    |
| aincs       | Array of angle increments.        |
| chords      | Array of chord lengths.           |
| xyzles      | Array of leading edge coordinates.|

<!-- ## Parameter sweep example


Initialize the AVLSolver:

```python
avl_solver = AVLSolver(geo_file="aircraft_mod.avl")
```
Retrieve the current geometry parameters:

```python

data = avl_solver.get_surface_params(
    include_geom=True,
    include_panneling=True,
    include_con_surf=True
)
```
Modify the desired parameters, e.g., scale the wing by a factor of 1.2:

```python

data["Wing"]["scale"] = np.array([1.2, 1.2, 1.2])
```
Set the modified parameters back to the AVL model:

```python

avl_solver.set_surface_params(data)
```
Execute the analysis run:

```python

avl_solver.set_constraint("alpha", 6.00)
avl_solver.execute_run()
```

Fetch the results and proceed with the next step in your parameter sweep. -->
