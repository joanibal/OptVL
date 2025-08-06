"""This script is intended to demonstrate how to use a custom geometry component with optvl"""
import openmdao.api as om
from optvl import OVLSolver, OVLGroup, OVLMeshReader
import numpy as np
import copy

class GeometryParametrizationComp(om.ExplicitComponent):
    def setup(self):
        # Input variables
        self.add_input('xles_in', shape_by_conn=True, desc='Baseline x leading edge coordinates')
        self.add_input('yles_in', shape_by_conn=True, desc='Baseline y leading edge coordinates')
        self.add_input('zles_in', shape_by_conn=True, desc='Baseline z leading edge coordinates')
        self.add_input('added_sweep', val=0.0, desc='added Sweep angle in degrees', units='deg')
        self.add_input('added_dihedral', val=0.0, desc='added Dihedral angle in degrees', units='deg')

        # Output variables
        self.add_output('xles_out', shape_by_conn=True, copy_shape='xles_in', desc='Transformed x leading edge coordinates')
        self.add_output('yles_out', shape_by_conn=True, copy_shape='yles_in', desc='Transformed y leading edge coordinates')
        self.add_output('zles_out', shape_by_conn=True, copy_shape='zles_in', desc='Transformed z leading edge coordinates')

        # Finite difference partials
        self.declare_partials("*", "*", method="cs")

    def compute(self, inputs, outputs):
        # Extracting input values
        xles_baseline = inputs['xles_in']
        yles_baseline = inputs['yles_in']
        zles_baseline = inputs['zles_in']
        transformed_xles = copy.deepcopy(xles_baseline)
        transformed_zles = copy.deepcopy(zles_baseline)
        dsweep = inputs['added_sweep']
        ddihedral = inputs['added_dihedral']
        
        dy = yles_baseline[-1] - yles_baseline[0]
        dx = xles_baseline[-1] - xles_baseline[0]
        
        sweep_baseline = np.arctan(dx/dy) 
        dx_new =  np.tan(sweep_baseline + np.pi/180 *dsweep) * dy
        ddx = dx_new - dx
        
        # linearly apply the change to the whole wing 
        transformed_xles +=  (yles_baseline - yles_baseline[0])/dy * ddx
        
        
        dz = zles_baseline[-1] - zles_baseline[0]
        
        dihedral_baseline = np.arctan(dz/dy) 
        dz_new =  np.tan(dihedral_baseline + np.pi/180 *ddihedral) * dy
        ddz = dz_new - dz
        
        # linearly apply the change to the whole wing 
        transformed_zles +=  (yles_baseline - yles_baseline[0])/dy * ddz
        
        outputs['xles_out'] = transformed_xles 
        outputs['yles_out'] = inputs['yles_in'] 
        outputs['zles_out'] = transformed_zles 
        


model = om.Group()
model.add_subsystem("mesh", OVLMeshReader(geom_file="aircraft.avl"))
model.add_subsystem('wing_param', GeometryParametrizationComp())
model.connect("mesh.Wing:xles",['wing_param.xles_in'] )
model.connect("mesh.Wing:yles",['wing_param.yles_in'] )
model.connect("mesh.Wing:zles",['wing_param.zles_in'] )
model.connect("wing_param.xles_out",['ovlsolver.Wing:xles'] )
model.connect("wing_param.yles_out",['ovlsolver.Wing:yles'] )
model.connect("wing_param.zles_out",['ovlsolver.Wing:zles'] )

model.add_subsystem("ovlsolver", OVLGroup(geom_file="aircraft.avl"))
model.add_design_var("ovlsolver.Wing:aincs", lower=-10, upper=10)
model.add_design_var("wing_param.added_sweep", lower=-10, upper=10)

# the outputs of AVL can be used as contraints
model.add_constraint("ovlsolver.CL", equals=1.5)
model.add_constraint("ovlsolver.Cm", equals=0.0, scaler=1e3)
# Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
# To help fix this add a contraint that keeps the variable monotonic

model.add_objective("ovlsolver.CD", scaler=1e2)

prob = om.Problem(model)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
prob.driver.options['tol'] = 1e-6
prob.driver.options['disp'] = True

prob.setup(mode='rev')
prob.run_driver()
om.n2(prob, show_browser=False, outfile="vlm_opt_param.html")


# do this instead if you want to check derivatives
# prob.run_model()
# prob.check_totals()


prob.model.ovlsolver.solver.avl.write_geom_file('opt_airplane.avl')
