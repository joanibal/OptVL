# Running alpha and CL sweeps

After initializing and setting up your `OVLSolver`, you can perform analysis tasks such as alpha and CL sweeps.
This is done by specifying new parameters and re-executing the run. 
The iterations over Mach number take longer to run, because the matrix defining the influence of the vortices needs to be reconstructed for each new Mach number. 

```python 
{% include "../examples/run_aero_sweeps.py" %}
```