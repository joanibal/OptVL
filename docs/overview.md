# Basic Usage Overview

`OptVL` offers a python API mirroring AVL's text interface.
The typical workflow involves loading a geometry file, adding constraints, and executing analysis runs.
In addition to this though, OptVL also offers unique optimization functionality. 

Below is an example of a basic analysis run to give you an idea of the workflow. 

```python 
{% include "../examples/run_basic.py" %}
```

The script will plot the geometry, set parameters, run an analysis and then show a CP surface plot.

The geometry plot should looke like this 
![aircraft geometry plot](figures/aircraft_geom.png)
and the CP surface plot should look like this
![aircraft cp](figures/aircraft_cp.png)
You can rotate, zoom, and pan in the window to look at different parts for the aircraft.
![aircraft cp view 2](figures/aircraft_cp_view2.png)

# Limitations
1. OptVL does not support multiple run cases
<!-- 2. There is no single precision version (8 digit output instead of 16 digits) of OptVL available for download
   - You could compile one yourself if you really need this -->
2. No support for the DESIGN keyword for setting twist in the avl geometry file since we use a different system for modifying all aspects of the geometry.
3. The windows can only use a maximum of 5400 vortices compared to the default of 6000. 
