# Optimization Overview

In general, you could connect OptVL to any gradient-free or gradient-based solver. 
In these docs, I focus on gradient-based method due to their greater efficiency. 
OptVL provides the derivatives, but to perform gradient-based optimization you need to pass those derivatives to an optimizer. 
However, here I'm going to show you the two ways I recommend connecting OptVL to an optimizer, [Scipy](https://scipy.org/) and [OpenMDAO](https://openmdao.org/). 

Between the two I recommend Scipy if you are doing something with few design variables and do not already know OpenMDAO. 
However, once you start incorporating more design variables and constraints or want to use a custom geometric parameterization you are probably better off learning OpenMDAO. 
If you are on the fence I recommend giving OpenMDAO a try. 
Your optimization problem will likely only grow more complex with time and OpenMDAO will make it easier to accommodate the growth in complexity.


See optimization setup with [Scipy see this page](optimization_setup_scipy.md) and to see optimization setup with [OpenMDAO see this page](optimization_setup_om.md)


## Tips
- Don't set trim constraints (e.g. CL = 0.5) on the OptVL solver instead set it as a constraint at the optimization level. The examples all show this. 
-  Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
   To help fix this, add a constraint that keeps the variable monotonic or use a custom parameterization.
- Discontinuities can appear when moving flaps or ailerons due to sparse paneling. Use section paneling for this case to preserve good paneling at the edges of the control surfaces.  