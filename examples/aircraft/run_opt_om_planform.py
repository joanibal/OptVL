"""A openmdao based optimization for an aicraft using optvl"""
import openmdao.api as om
from optvl import AVLGroup, Differencer, AVLMeshReader
import numpy as np
import copy


class GeometryParametrizationComp(om.ExplicitComponent):
    def setup(self):
        # Input variables
        self.add_input('xles_in', shape_by_conn=True, desc='Baseline x leading edge coordinates')
        self.add_input('yles_in', shape_by_conn=True, desc='Baseline y leading edge coordinates')
        self.add_input('added_sweep', val=0.0, desc='added Sweep angle in degrees', units='deg')

        # Output variables
        self.add_output('xles_out', shape_by_conn=True, copy_shape='xles_in', desc='Transformed xyz leading edge coordinates')

        # Finite difference partials
        self.declare_partials("*", "*", method="cs")

    def compute(self, inputs, outputs):
        # Extracting input values
        yles = inputs['yles_in']
        xles = inputs['xles_in']
        span = yles[-1]
        relative_span = yles/span
        
        outputs['xles_out'] = xles + inputs['added_sweep']*relative_span


model = om.Group()
geom_dvs = model.add_subsystem("geom_dvs", om.IndepVarComp())

geom_dvs.add_output('chords', shape_by_conn=True)
model.connect('geom_dvs.chords', "avlsolver.Wing:chords")
geom_dvs.add_output('aincs', shape_by_conn=True)
model.connect('geom_dvs.aincs', "avlsolver.Wing:aincs")


model.add_subsystem("mesh", AVLMeshReader(geom_file="rectangle.avl"))
model.add_subsystem('geom_param', GeometryParametrizationComp())

model.connect("mesh.Wing:xles",['geom_param.xles_in'] )
model.connect("geom_param.xles_out",['avlsolver.Wing:xles'] )

model.connect("mesh.Wing:yles",['geom_param.yles_in'] )

model.add_subsystem("avlsolver", AVLGroup(geom_file="rectangle.avl", output_stabililty_derivs=True, write_grid=True, output_dir='opt_output_sweep'))


model.add_subsystem("differ_chords", Differencer())
model.connect("geom_dvs.chords", "differ_chords.input_vec")

model.add_subsystem("differ_aincs", Differencer())
model.connect("geom_dvs.aincs", "differ_aincs.input_vec")

# look at vlm_opt.html to see all the design variables and add them here
model.add_design_var("geom_param.added_sweep", lower=0.0, upper=0.66)

model.add_design_var("avlsolver.Wing:aincs", lower=-10, upper=10)
# model.add_design_var("avlsolver.Elevator", lower=-10, upper=10)
# model.add_design_var("avlsolver.Wing:scale", indices=[1])
model.add_design_var("avlsolver.Wing:chords", lower=0.1, upper=1.0)


# the outputs of AVL can be used as contraints
model.add_constraint("avlsolver.CL", equals=1.5)
model.add_constraint("avlsolver.CM", equals=0.0, scaler=1e3)
model.add_constraint("avlsolver.dCM/dalpha", lower=-0.2, upper=0.0, scaler=1e3)


# Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
# To help fix this add a contraint that keeps the variable monotonic
model.add_constraint("differ_aincs.diff_vec",  upper=0.0, linear=True)
model.add_constraint("differ_chords.diff_vec", lower=-0.33, upper=0.0, linear=True)

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
# prob.run_model()

# aincs = prob.get_val('avlsolver.Wing:aincs')
# print(f'avlsolver.Wing:aincs {aincs}')
# del_ele = prob.get_val('avlsolver.Elevator')
# print(f'avlsolver.Elevator {del_ele}')
# do this instead if you want to check derivatives
# prob.run_model()
# prob.check_totals()
