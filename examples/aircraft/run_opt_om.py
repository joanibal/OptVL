"""A openmdao based optimization for an aicraft using optvl"""
import openmdao.api as om
from optvl import AVLGroup


model = om.Group()
model.add_subsystem("avlsolver", AVLGroup(geom_file="aircraft.avl"))
# look at vlm_opt.html to see all the design variables and add them here
model.add_design_var("avlsolver.Wing:aincs", lower=-10, upper=10)
model.add_design_var("avlsolver.Elevator", lower=-10, upper=10)

# the outputs of AVL can be used as contraints
model.set_constraint("avlsolver.CL", equals=1.5)
model.set_constraint("avlsolver.CM", equals=0.0, scaler=1e3)
# Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
# To help fix this add a contraint that keeps the variable monotonic

# the scaler values bring the objective functinon to ~ order 1 for the optimizer
model.add_objective("avlsolver.CD", scaler=1e3)

prob = om.Problem(model)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
prob.driver.options['tol'] = 1e-6
prob.driver.options['disp'] = True

prob.setup(mode='rev')
om.n2(prob, show_browser=False, outfile="vlm_opt.html")
prob.run_driver()

aincs = prob.get_val('avlsolver.Wing:aincs')
print(f'avlsolver.Wing:aincs {aincs}')
del_ele = prob.get_val('avlsolver.Elevator')
print(f'avlsolver.Elevator {del_ele}')
# do this instead if you want to check derivatives
# prob.run_model()
# prob.check_totals()
