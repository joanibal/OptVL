# Visualization


## With Matplotlib/python
The visualizations described in this section rely on the `matplotlib` python package, but no external tools. 

### Looking at the geometry
To quickly look at your geometry you can use the
```python
ovl.plot_geom()
``` 
command, which will produce a figure that looks like this. 
![aircraft geometry plot](figures/aircraft_geom.png)


### Cp plots 
To get a quick view of the coefficient of pressure on the surfaces of the aircraft you can use 
```python
ovl.plot_cp()
```
Which will produce a plot like this. 

![aircraft cp](figures/aircraft_cp.png)
You can rotate, zoom, and pan in the window to look at different parts for the aircraft.
![aircraft cp view 2](figures/aircraft_cp_view2.png)

## Paraview or tecplot
This section describes how to produce data that ca be used with either tecplot or Paraview.

To generate a tecplot ascii file (which can be read by both tecplot and paraview). 
After and analysis use the command 
```python
ovl_solver.write_tecplot('test')
```
. This command will write a file called `test.dat`.

### Tecplot 360

If you load the file with Tecplot 360 it will automatically be detected as a Tecplot acii file and will be loaded accordingly. 
By default Tecplot tends to plot the data with a 2D Cartesian plot, but you will have to switch it 3D Cartesian to see the full geometry. 

![](figures/tecplot_geom.png)

You can also plot the mesh and the coefficient of pressure on the surface. 

![](figures/tecplot_cp.png)

When writing out the tecplot file you can also specify a solution time with the optional keyword argument. 
This makes it easier to flip through different data files in tecplot but makes the output incompatible with paraview.
```python 
ovl_solver.write_tecplot('test', solution_time=1)
```

## Paraview

When loading your data into paraview be sure to use the tecplot dat loader. 

![](figures/paraview_loader.png)

With the data loaded you can view the 3D shape, Coefficient of Pressure, and mesh just like with Tecplot. 

![](figures/paraview_cp.png)


Paraview nicely stacks up files of the same name but with different numeric suffixes and loads them into different timesteps. 
So if you want to easily view the results of an optimization or parameter sweep name your output files the same name with a different number for the suffix. 