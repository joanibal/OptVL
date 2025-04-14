# Optimization 

See [`examples/aircraft/run_opt.py`](https://github.com/joanibal/OptAVL/blob/master/examples/aircraft/run_opt.py) for a basic example. 
Other scripts in the examples' directory provide examples for more complex cases. 

More documentation to follow. 
For now if you have a question post an issue on github.  


## Debugging
-  Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
   To help fix this, add a constraint that keeps the variable monotonic or use a custom parameterization.
- Discontinuities can appear when moving flaps or ailerons due to sparse paneling. Use section paneling for this case to preserve good paneling at the edges of the control surfaces.  