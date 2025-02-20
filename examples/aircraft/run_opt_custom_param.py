"""This script is intended to demonstrate how to use a custom geometry component with optvl"""
import openmdao.api as om
from optvl import AVLSolver, AVLGroup, AVLMeshReader
import numpy as np
import copy

class GeometryParametrizationComp(om.ExplicitComponent):
    def setup(self):
        # Input variables
        self.add_input('xyzles_in', shape_by_conn=True, desc='Baseline xyz leading edge coordinates')
        self.add_input('added_sweep', val=0.0, desc='added Sweep angle in degrees', units='deg')
        self.add_input('added_dihedral', val=0.0, desc='added Dihedral angle in degrees', units='deg')

        # Output variables
        self.add_output('xyzles_out', shape_by_conn=True, copy_shape='xyzles_in', desc='Transformed xyz leading edge coordinates')

        # Finite difference partials
        self.declare_partials("*", "*", method="cs")

    def compute(self, inputs, outputs):
        # Extracting input values
        xyzles_baseline = inputs['xyzles_in']
        transformed_xyzles = copy.deepcopy(xyzles_baseline)
        dsweep = inputs['added_sweep']
        ddihedral = inputs['added_dihedral']
        
        dy = xyzles_baseline[-1, 1] - xyzles_baseline[0, 1]
        dx = xyzles_baseline[-1, 0] - xyzles_baseline[0, 0]
        
        sweep_baseline = np.arctan(dx/dy) 
        dx_new =  np.tan(sweep_baseline + np.pi/180 *dsweep) * dy
        ddx = dx_new - dx
        
        # linearly apply the change to the whole wing 
        transformed_xyzles[:, 0] +=  (xyzles_baseline[:, 1] - xyzles_baseline[0, 1])/dy * ddx
        
        
        dz = xyzles_baseline[-1, 2] - xyzles_baseline[0, 2]
        
        dihedral_baseline = np.arctan(dz/dy) 
        dz_new =  np.tan(dihedral_baseline + np.pi/180 *ddihedral) * dy
        ddz = dz_new - dz
        
        # linearly apply the change to the whole wing 
        transformed_xyzles[:, 2] +=  (xyzles_baseline[:, 1] - xyzles_baseline[0, 1])/dy * ddz
        
        outputs['xyzles_out'] = transformed_xyzles


model = om.Group()
model.add_subsystem("mesh", AVLMeshReader(geom_file="aircraft.avl"))
model.add_subsystem('wing_param', GeometryParametrizationComp())
model.connect("mesh.Wing:xyzles",['wing_param.xyzles_in'] )
model.connect("wing_param.xyzles_out",['avlsolver.Wing:xyzles'] )

model.add_subsystem("avlsolver", AVLGroup(geom_file="aircraft.avl"))
model.add_design_var("avlsolver.Wing:aincs", lower=-10, upper=10)
model.add_design_var("wing_param.added_sweep", lower=-10, upper=10)

# the outputs of AVL can be used as contraints
model.set_constraint("avlsolver.CL", equals=1.5)
model.set_constraint("avlsolver.CM", equals=0.0, scaler=1e3)
# Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
# To help fix this add a contraint that keeps the variable monotonic

model.add_objective("avlsolver.CD", scaler=1e2)

prob = om.Problem(model)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
prob.driver.options['tol'] = 1e-6
prob.driver.options['disp'] = True

prob.setup(mode='rev')
om.n2(prob, show_browser=False, outfile="vlm_opt_param.html")
prob.run_driver()

# do this instead if you want to check derivatives
# prob.run_model()
# prob.check_totals()


prob.model.avlsolver.solver.avl.write_geom_file('opt_airplane.avl')
