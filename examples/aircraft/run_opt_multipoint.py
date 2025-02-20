"""A openmdao based optimization for an aicraft using optvl"""
import openmdao.api as om
from optvl import AVLSolver, AVLGroup
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--record", default=False, action="store_true")
parser.add_argument("--plot_opt_hist", default=False, action="store_true")
args = parser.parse_args()

class Top(om.Group):

    def setup(self):    
        # add independent var comps to hold design variables
        self.add_subsystem("flt_cond_dvs", om.IndepVarComp())
        self.add_subsystem("geom_dvs", om.IndepVarComp())
        self.add_subsystem("con_surf_dvs", om.IndepVarComp())
        
        self.add_subsystem("perturbed_flt_cond", om.ExecComp(['ptb_alpha=alpha+0.1'], units='deg'))
        
        self.add_subsystem("cruise_avl", AVLGroup(geom_file="aircraft.avl", mass_file="aircraft.mass", output_stabililty_derivs=True))
        self.add_subsystem("perturbed_cruise_avl", AVLGroup(geom_file="aircraft.avl", mass_file="aircraft.mass"))
        
    def configure(self):
        # add the flight flt_condition inputs as possible design variables
        # BUT perturb the alpha passed into perturbed_cruise_avl
        self.flt_cond_dvs.add_output('alpha', val=0.0, units='deg')
        self.flt_cond_dvs.add_output('beta', val=0.0, units='deg')
        self.connect('flt_cond_dvs.alpha', 'cruise_avl.alpha')
        self.connect('flt_cond_dvs.beta', ['cruise_avl.beta', 'perturbed_cruise_avl.beta'])
        
        self.connect('flt_cond_dvs.alpha', 'perturbed_flt_cond.alpha')
        self.connect('perturbed_flt_cond.ptb_alpha', 'perturbed_cruise_avl.alpha')
        
        # add all the geometric inputs as possible design variables 
        # Connect the values of each cruise compoenent
        geom_inputs = self.cruise_avl.solver.list_inputs(tags="geom", val=True, units=True, out_stream=None)
        for geom_input in geom_inputs:
            meta_data = geom_input[1]
            units = meta_data['units']
            promoted_name = meta_data['prom_name']
            # get the initail values from the component otherwise it would be set to all 1's
            val = meta_data['val']
            self.geom_dvs.add_output(promoted_name, val=val, units=units)
            
            cruise_name = f"{self.cruise_avl.name}.{promoted_name}"
            ptb_cruise_name = f"{self.perturbed_cruise_avl.name}.{promoted_name}"
            self.connect(f'geom_dvs.{promoted_name}', cruise_name )
            self.connect(f'geom_dvs.{promoted_name}', ptb_cruise_name )
        
        # add all the control surface inputs as possible design variables 
        # Connect the values of each cruise compoenent
        con_surf_inputs = self.cruise_avl.solver.list_inputs(tags="con_surf", val=True, units=True, out_stream=None)
        for con_surf_input in con_surf_inputs:
            meta_data = con_surf_input[1]
            units = meta_data['units']
            promoted_name = meta_data['prom_name']
            # get the initail values from the component otherwise it would be set to all 1's
            val = meta_data['val']
            self.con_surf_dvs.add_output(promoted_name, val=val, units=units)

            cruise_name = f"{self.cruise_avl.name}.{promoted_name}"
            ptb_cruise_name = f"{self.perturbed_cruise_avl.name}.{promoted_name}"
            self.connect(f'con_surf_dvs.{promoted_name}', cruise_name )
            self.connect(f'con_surf_dvs.{promoted_name}', ptb_cruise_name)
        
model = Top()

# look at vlm_mp_opt.html to see all the design variables and add them here
# tip dihedral  
idx_sec = 4  # idx of the last section
idx_xyz = 2  # z values are in the third place
flat_idx = 3*idx_sec + idx_xyz
model.add_design_var("geom_dvs.Wing:xyzles", indices=[flat_idx], flat_indices=True,  lower=0, upper=1)
model.add_design_var("geom_dvs.Wing:aincs", lower=-10, upper=10)
model.add_design_var("flt_cond_dvs.alpha", lower=-10, upper=10)
model.add_design_var("con_surf_dvs.Elevator", lower=-10, upper=10)

# the outputs of AVL can be used as contraints
model.set_constraint("cruise_avl.CL", equals=1.5)
model.set_constraint("cruise_avl.CM", equals=0.0)

model.set_constraint("perturbed_cruise_avl.CM", upper=-1e-3) # make sure that dCM_dAlpha is less than 0 for static stability

model.add_objective("cruise_avl.CD", ref=1e-2)
# Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
# To help fix this add a contraint that keeps the variable monotonic

prob = om.Problem(model)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['debug_print'] = ["desvars", "ln_cons", "nl_cons", "objs"]
prob.driver.options['tol'] = 1e-6
prob.driver.options['disp'] = True
prob.driver.options['maxiter'] = 100

if args.record:
    case_recorder_file = "om_opt_hist.db"
    recorder = om.SqliteRecorder(case_recorder_file)
    prob.driver.add_recorder(recorder)

    prob.driver.recording_options["record_desvars"] = True
    prob.driver.recording_options["record_responses"] = True
    prob.driver.recording_options["record_objectives"] = True
    prob.driver.recording_options["record_constraints"] = True

prob.setup(mode='rev')
om.n2(prob, show_browser=False, outfile="vlm_mp_opt.html")

prob.run_driver()
# do this instead if you want to check derivatives
# prob.run_model()
# prob.check_totals()

prob.model.cruise_avl.solver.avl.write_geom_file('opt_airplane_mp.avl')

if args.plot_opt_hist and args.record:
    import matplotlib.pyplot as plt

    cr = om.CaseReader(case_recorder_file)

    driver_cases = cr.list_cases('driver', out_stream=None)

    num_iters = len(driver_cases)

    case = cr.get_case(driver_cases[0])

    obj_data = case.get_objectives()
    con_data = case.get_constraints()


    # set the contraint values to save
    case = cr.get_case(driver_cases[0])

    for idx_iter in range(1,num_iters):
        case = cr.get_case(driver_cases[idx_iter])
        obj_data_iter = case.get_objectives()
        for obj_key in obj_data_iter:
            obj_data[obj_key] = np.append(obj_data[obj_key], obj_data_iter[obj_key])
            print(idx_iter,obj_key, obj_data_iter[obj_key])
        
        con_data_iter = case.get_constraints()
        for con_key in con_data_iter:
            con_data[con_key] = np.append(con_data[con_key], con_data_iter[con_key])
            print(idx_iter,con_key, con_data_iter[con_key])

    # Prepare x-axis data - iterations
    iterations = np.arange(num_iters)

    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # Plotting objectives
    for key, values in obj_data.items():
        ax1.plot(iterations, values, label=key)

    # Adding legend and labels
    ax1.set_ylabel('Objective Values')
    ax1.legend()
    ax1.set_title('Objective Data Over Iterations')

    # Plotting constraints
    for key, values in con_data.items():
        ax2.plot(iterations, values, label=key)

    # Adding legend and labels
    ax2.set_ylabel('Constraint Values')
    ax2.legend()
    ax2.set_title('Constraint Data Over Iterations')

    # Setting common x-axis label
    plt.xlabel('Iterations')

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()

