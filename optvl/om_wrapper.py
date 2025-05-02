import os
import openmdao.api as om
from optvl import OVLSolver
import numpy as np
import copy
import time

class OVLGroup(om.Group):
    """This is the main the top level group for interacting with OptVL. 

    Args:
        geom_file: the input geometry file
        mass_file: the optional mass file  
        write_grid: should the grid be written after avery analysis 
        write_grid_sol_time: add the iteration count as the solution time for easier postprocessing in tecplot
        output_dir: the output directory for the files generated
        input_param_vals: flag to turn on the flght parameters (Mach, Velocity, etc.) as inputs
        input_ref_val: flag to turn on the geometric reference values (Sref, Cref, Bref) as inputs
        output_stabililty_derivs: flag to turn on the output of stability derivatives
        output_con_surf_derivs: flag to turn on the output of control surface deflections 
    """
    def initialize(self):
        self.options.declare("geom_file", types=str)
        self.options.declare("mass_file", default=None)
        self.options.declare("write_grid", types=bool, default=False)
        self.options.declare("write_grid_sol_time", types=bool, default=False)
        self.options.declare("output_dir", types=str, recordable=False, default=".")

        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        
        self.options.declare("output_stabililty_derivs", types=bool, default=False)
        self.options.declare("output_con_surf_derivs", types=bool, default=False)

    def setup(self):
        geom_file = self.options["geom_file"]
        mass_file = self.options["mass_file"]
        
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        
        output_stabililty_derivs = self.options["output_stabililty_derivs"]
        output_con_surf_derivs = self.options["output_con_surf_derivs"]

        ovl = OVLSolver(geo_file=geom_file, mass_file=mass_file, debug=False)

        self.add_subsystem("solver", OVLSolverComp(ovl=ovl, 
                                    input_param_vals=input_param_vals,
                                    input_ref_vals=input_ref_vals),
                                    promotes=["*"])
        self.add_subsystem("funcs", OVLFuncsComp(ovl=ovl,
                                    input_param_vals=input_param_vals,
                                    input_ref_vals=input_ref_vals,
                                    output_stabililty_derivs=output_stabililty_derivs,
                                    output_con_surf_derivs=output_con_surf_derivs),
                                    promotes=["*"])
        if self.options["write_grid"]:
            self.add_subsystem(
                "postprocess", OVLPostProcessComp(ovl=ovl, 
                                                    output_dir=self.options["output_dir"],
                                                    write_grid_sol_time=self.options["write_grid_sol_time"],
                                                    input_param_vals=input_param_vals,
                                                    input_ref_vals=input_ref_vals),
                                                    promotes=["*"]
                                                )

# helper functions used by the AVL components
def add_ovl_controls_as_inputs(self, ovl):
    # add the control surfaces as inputs
    self.control_names = ovl.get_control_names()
    for c_name in self.control_names:
        self.add_input(c_name, val=0.0, units="deg", tags="con_surf")
    return self.control_names

def add_ovl_geom_vars(self, ovl, add_as="inputs"):
    # add the geometric parameters as inputs
    surf_data = ovl.get_surface_params()

    for surf in surf_data:
        for key in surf_data[surf]:
            geom_key = f"{surf}:{key}"
            if add_as == "inputs":
                self.add_input(geom_key, val=surf_data[surf][key], tags="geom")
            elif add_as == "outputs":
                self.add_output(geom_key, val=surf_data[surf][key], tags="geom")
                
def add_ovl_conditions_as_inputs(sys, ovl):
    # TODO: add all the condition constraints
        
    sys.add_input("alpha", val=0.0, units="deg", tags="flt_cond")
    sys.add_input("beta", val=0.0, units="deg" , tags="flt_cond")


def add_ovl_params_as_inputs(sys, ovl):
    # TODO: add all par vals with the analysis is supported
    
    # only adding the ones people would use for now
    for param in ["CD0", "Mach", "X cg", "Y cg", "Z cg"]:
        val = ovl.get_parameter(param)
        sys.add_input(param, val=val, tags="param")

def add_ovl_refs_as_inputs(sys, ovl):
    ref_data = ovl.get_reference_data()
    
    for key, val in ref_data.items():
        sys.add_input(key, val=val, tags="ref_val")


def om_input_to_surf_dict(sys, inputs):
    geom_inputs = sys.list_inputs(tags="geom", val=False, out_stream=None)
    # convert to a list witout tuples
    geom_inputs = [x[0] for x in geom_inputs]

    surf_data = {}
    for input_var in inputs:
        if input_var in geom_inputs:
            # split the input name into surface name and parameter name
            surf, param = input_var.split(":")

            # update the corresponding parameter in the surface data
            if surf not in surf_data:
                surf_data[surf] = {}
            if inputs[input_var].size == 1:
                # if the input is a scalar, convert to a 1D array for numpy
                surf_data[surf][param] = inputs[input_var][0]
            else:
                surf_data[surf][param] = inputs[input_var]
    return surf_data


def om_surf_dict_to_input(surf_dict):
    input_data = {}
    for surf_key in surf_dict:
        for geom_key in surf_dict[surf_key]:
            input_var = ":".join([surf_key, geom_key])

            input_data[input_var] = surf_dict[surf_key][geom_key]

    return input_data

def om_set_avl_inputs(sys, inputs):
    for c_name in sys.control_names:
        sys.ovl.set_constraint(c_name, inputs[c_name][0])

    sys.ovl.set_constraint("alpha", inputs["alpha"][0])
    sys.ovl.set_constraint("beta", inputs["beta"][0])

    # add the parameters to the run
    for param in sys.ovl.param_idx_dict:
        if param in inputs:
            val = inputs[param][0]
            sys.ovl.set_parameter(param, val)
    
    # add the parameters to the run
    for ref in sys.ovl.ref_var_to_fort_var:
        if ref in inputs:
            val = inputs[ref][0]
            sys.ovl.set_reference_data({ref: val})
    


class OVLSolverComp(om.ImplicitComponent):
    """
    OpenMDAO component that wraps optvl solver. This is added as part of the OVLgroup
    """

    def initialize(self):
        self.options.declare("ovl", types=OVLSolver, recordable=False)
        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        

    def setup(self):
        self.ovl = self.options["ovl"]
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        
        self.num_states = self.ovl.get_mesh_size()
        self.num_cs = self.ovl.get_num_control_surfs()
        self.num_vel = self.ovl.NUMAX

        self.add_output("gamma", val=np.zeros(self.num_states))
        self.add_output("gamma_d", val=np.zeros((self.num_cs, self.num_states)))
        self.add_output("gamma_u", val=np.zeros((self.num_vel, self.num_states)))
        
        add_ovl_conditions_as_inputs(self, self.ovl)
        
        if input_param_vals:
            add_ovl_params_as_inputs(self, self.ovl)
        
        if input_ref_vals:
            add_ovl_refs_as_inputs(self, self.ovl)
        
        self.control_names = add_ovl_controls_as_inputs(self, self.ovl)
        add_ovl_geom_vars(self, self.ovl, add_as="inputs")
        
        self.res_slice = (slice(0, self.num_states),)
        self.res_d_slice = (slice(0, self.num_cs), slice(0, self.num_states))
        self.res_u_slice = (slice(0, self.num_vel), slice(0, self.num_states))


    def apply_nonlinear(self, inputs, outputs, residuals):
        om_set_avl_inputs(self, inputs)
        
        surf_data = om_input_to_surf_dict(self, inputs)
        self.ovl.set_surface_params(surf_data)

        gam_arr = outputs["gamma"]
        gam_d_arr = outputs["gamma_d"]
        gam_u_arr = outputs["gamma_u"]


        
        # TODO-api: this should probably be an API level call to set gamma
        self.ovl.set_avl_fort_arr("VRTX_R", "GAM", gam_arr, slicer=self.res_slice)
        self.ovl.set_avl_fort_arr("VRTX_R", "GAM_D", gam_d_arr, slicer=self.res_d_slice)
        self.set_avl_fort_arr("VRTX_R", "GAM_U", gam_u_arr, slicer=self.res_u_slice)

        # propogate the seeds through without resolving
        self.ovl.avl.update_surfaces()
        self.ovl.avl.get_res()


        res = copy.deepcopy(self.ovl.get_avl_fort_arr("VRTX_R", "RES", slicer=self.res_slice))
        residuals["gamma"] = res
        res_d = copy.deepcopy(self.ovl.get_avl_fort_arr("VRTX_R", "RES_D", slicer=self.res_d_slice))
        residuals["gamma_d"] = res_d
        res_u = copy.deepcopy(self.ovl.get_avl_fort_arr("VRTX_R", "RES_U", slicer=self.res_u_slice))
        residuals["gamma_u"] = res_u
        
        # this routine shouldn't be used normally

    def solve_nonlinear(self, inputs, outputs):
        start_time = time.time()
        om_set_avl_inputs(self, inputs)

        # update the surface parameters
        surf_data = om_input_to_surf_dict(self, inputs)
        self.ovl.set_surface_params(surf_data)

        # def_dict = self.ovl.get_control_deflections()
        print('executing ovl run')
        self.ovl.execute_run()

        gam_arr = self.ovl.get_avl_fort_arr("VRTX_R", "GAM", slicer=self.res_slice)

        outputs["gamma"] = copy.deepcopy(gam_arr)

        gam_d_arr = self.ovl.get_avl_fort_arr("VRTX_R", "GAM_D", slicer=self.res_d_slice)
        outputs["gamma_d"] = copy.deepcopy(gam_d_arr)
        
        gam_u_arr = self.ovl.get_avl_fort_arr("VRTX_R", "GAM_U", slicer=self.res_u_slice)
        outputs["gamma_u"] = copy.deepcopy(gam_u_arr)

        # run_data = self.ovl.get_total_forces()
        # for func_key in run_data:
        #     print(func_key, run_data[func_key])
        # func_key = "CL"
        # print(func_key, run_data[func_key])
        # print("AVL solve time: ", time.time() - start_time)

    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):
        if mode == "fwd":
            con_seeds = {}
            for con_key in ["alpha", "beta"]:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
            for con_key in self.control_names:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
                    
            geom_seeds = om_input_to_surf_dict(self, d_inputs)
            
            
            param_seeds = {}
            for param in self.ovl.param_idx_dict:
                if param in d_inputs:
                    param_seeds[param] = d_inputs[param]
                
            ref_seeds = {}
            for ref in self.ovl.ref_var_to_fort_var:
                if ref in d_inputs:
                    ref_seeds[ref] = d_inputs[ref]

            _, res_seeds, _, _, res_d_seeds, res_u_seeds = self.ovl._execute_jac_vec_prod_fwd(con_seeds=con_seeds, geom_seeds=geom_seeds,
                                                                           param_seeds=param_seeds, ref_seeds=ref_seeds)

            d_residuals["gamma"] += res_seeds
            d_residuals["gamma_d"] += res_d_seeds
            d_residuals["gamma_u"] += res_u_seeds

        if mode == "rev":
            if "gamma" in d_residuals:
                self.ovl.clear_ad_seeds_fast()
                res_seeds = d_residuals["gamma"]
                res_d_seeds = d_residuals["gamma_d"]
                res_u_seeds = d_residuals["gamma_u"]

                con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds, gamma_u_seeds, param_seeds, ref_seeds = self.ovl._execute_jac_vec_prod_rev(
                    res_seeds=res_seeds, res_d_seeds=res_d_seeds, res_u_seeds=res_u_seeds
                )

                if "gamma" in d_outputs:
                    d_outputs["gamma"] += gamma_seeds

                if "gamma_d" in d_outputs:
                    d_outputs["gamma_d"] += gamma_d_seeds
                
                if "gamma_u" in d_outputs:
                    d_outputs["gamma_u"] += gamma_u_seeds

                d_input_geom = om_surf_dict_to_input(geom_seeds)

                for d_input in d_inputs:
                    if d_input in d_input_geom:
                        d_inputs[d_input] += d_input_geom[d_input]
                    elif d_input in ["alpha", "beta"]:
                        d_inputs[d_input] += con_seeds[d_input]
                    elif d_input in self.control_names:
                        d_inputs[d_input] += con_seeds[d_input]
                    elif d_input in param_seeds:
                        d_inputs[d_input] += param_seeds[d_input]
                    elif d_input in ref_seeds:
                        d_inputs[d_input] += ref_seeds[d_input]

                  

    def solve_linear(self, d_outputs, d_residuals, mode):
        if mode == "rev":
            self.ovl.set_gamma_ad_seeds(d_outputs["gamma"])
            self.ovl.set_gamma_d_ad_seeds(d_outputs["gamma_d"])
            self.ovl.set_gamma_u_ad_seeds(d_outputs["gamma_u"])
            # start_time = time.time()
            solve_stab_deriv_adj=True
            solve_con_surf_adj=True
            self.ovl.avl.solve_adjoint(solve_stab_deriv_adj, solve_con_surf_adj)
            # print("OM Solve adjoint time: ", time.time() - start_time)
            d_residuals["gamma"] = self.ovl.get_residual_ad_seeds()
            d_residuals["gamma_d"] = self.ovl.get_residual_d_ad_seeds()
            d_residuals["gamma_u"] = self.ovl.get_residual_u_ad_seeds()

        elif mode == "fwd":
            raise NotImplementedError("only reverse mode derivaties implemented. Use prob.setup(mode='rev')")


class OVLFuncsComp(om.ExplicitComponent):
    """This component uses OptVL to compute functionals given a circulation solution. This is added as part of the OVL group

    Args:
        om: _description_
    """
    def initialize(self):
        self.options.declare("ovl", types=OVLSolver, recordable=False)
        self.options.declare("output_stabililty_derivs", types=bool, default=False)
        self.options.declare("output_con_surf_derivs", types=bool, default=False)
        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        

    def setup(self):
        self.ovl = self.options["ovl"]
        self.num_states = self.ovl.get_mesh_size()
        self.num_cs = self.ovl.get_num_control_surfs()
        self.num_vel = self.ovl.NUMAX
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        
        self.add_input("gamma", val=np.zeros(self.num_states))
        self.add_input("gamma_d", val=np.zeros((self.num_cs, self.num_states)))
        self.add_input("gamma_u", val=np.zeros((self.num_vel, self.num_states)))

        add_ovl_conditions_as_inputs(self, self.ovl)

        if input_param_vals:
            add_ovl_params_as_inputs(self, self.ovl)
        
        if input_ref_vals:
            add_ovl_refs_as_inputs(self, self.ovl)
            
        self.control_names = add_ovl_controls_as_inputs(self, self.ovl)
        add_ovl_geom_vars(self, self.ovl, add_as="inputs")

        # add the outputs
        for func_key in self.ovl.case_var_to_fort_var:
            self.add_output(func_key)

        self.output_con_surf_derivs = self.options["output_con_surf_derivs"]
        if self.output_con_surf_derivs:
            for func_key in self.ovl.case_derivs_to_fort_var:
                self.add_output(func_key)
        
        self.output_stabililty_derivs = self.options["output_stabililty_derivs"]
        if self.output_stabililty_derivs:
            deriv_dict = self.ovl.case_stab_derivs_to_fort_var
            for func_key in deriv_dict:
                self.add_output(func_key)        
        
        # TODO-refactor: push these slices down into the ovl class?                     
        self.res_slice = (slice(0, self.num_states),)
        self.res_d_slice = (slice(0, self.num_cs), slice(0, self.num_states))
        self.res_u_slice = (slice(0, self.num_vel), slice(0, self.num_states))

    def compute(self, inputs, outputs):
        # self.ovl.set_gamma(inputs['gamma'])
        # for c_name in self.control_names:
        #     self.ovl.set_constraint(c_name, inputs[c_name][0])
        # def_dict = self.ovl.get_control_deflections()

        # TODO: set_constraint does not correctly do derives yet
        start_time = time.time()
        om_set_avl_inputs(self, inputs)

        # update the surface parameters
        surf_data = om_input_to_surf_dict(self, inputs)
        self.ovl.set_surface_params(surf_data)

        gam_arr = inputs["gamma"]
        gam_d_arr = inputs["gamma_d"]
        gam_u_arr = inputs["gamma_u"]

        self.ovl.set_avl_fort_arr("VRTX_R", "GAM", gam_arr, slicer=self.res_slice)
        self.ovl.set_avl_fort_arr("VRTX_R", "GAM_D", gam_d_arr, slicer=self.res_d_slice)
        self.ovl.set_avl_fort_arr("VRTX_R", "GAM_U", gam_u_arr, slicer=self.res_u_slice)
        
        # TODO: only update what you need to. 
        # residuals (and AIC?) do not need to be calculated 
        # but in get_res, alpha and beta are set
        self.ovl.avl.update_surfaces()
        self.ovl.avl.get_res()
        self.ovl.avl.velsum()
        self.ovl.avl.aero()

        run_data = self.ovl.get_total_forces()

        for func_key in run_data:
            # print(f' {func_key} {run_data[func_key]}')
            outputs[func_key] = run_data[func_key]
        # print(f" CD {run_data['CD']} CL {run_data['CL']}")

        if self.output_con_surf_derivs:
            consurf_derivs_seeds = self.ovl.get_control_stab_derivs()
            for func_key in consurf_derivs_seeds:
                    # var_name = f"d{func_key}_d{con_name}"
                outputs[func_key] = consurf_derivs_seeds[func_key]

        if self.output_stabililty_derivs:
            stab_derivs = self.ovl.get_stab_derivs()
            for func_key in stab_derivs:
                outputs[func_key] = stab_derivs[func_key]
                
        # print("Funcs Compute time: ", time.time() - start_time)

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            con_seeds = {}
            for con_key in ["alpha", "beta"]:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
            
            for con_key in self.control_names:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
                    
            if "gamma" in d_inputs:
                gamma_seeds = d_inputs["gamma"]
            else:
                gamma_seeds = None

            if "gamma_d" in d_inputs:
                gamma_d_seeds = d_inputs["gamma_d"]
            else:
                gamma_d_seeds = None
                
            if "gamma_u" in d_inputs:
                gamma_u_seeds = d_inputs["gamma_u"]
            else:
                gamma_u_seeds = None
                
            param_seeds = {}
            for param in self.ovl.param_idx_dict:
                if param in d_inputs:
                    param_seeds[param] = d_inputs[param]
                
            ref_seeds = {}
            for ref in self.ovl.ref_var_to_fort_var:
                if ref in d_inputs:
                    ref_seeds[ref] = d_inputs[ref]

            geom_seeds = self.om_input_to_surf_dict(self, d_inputs)

            func_seeds, _, csd_seeds, stab_derivs_seeds, _, _ = self.ovl._execute_jac_vec_prod_fwd(
                con_seeds=con_seeds, geom_seeds=geom_seeds, gamma_seeds=gamma_seeds, gamma_d_seeds=gamma_d_seeds, gamma_u_seeds=gamma_u_seeds,
                param_seeds=param_seeds, ref_seeds=ref_seeds
            )

            for func_key in func_seeds:
                d_outputs[func_key] += func_seeds[func_key]

            for func_key in csd_seeds:
                for con_name in csd_seeds[func_key]:
                    var_name = f"d{func_key}_d{con_name}"
                    d_outputs[var_name] = csd_seeds[func_key][con_name]
            
            for func_key in stab_derivs_seeds:
                for var in stab_derivs_seeds[func_key]:
                    var_name = f"d{func_key}_d{var}"
                    d_outputs[var_name] = stab_derivs_seeds[func_key][var]

        if mode == "rev":
            self.ovl.clear_ad_seeds_fast()

            # add the outputs
            func_seeds = {}
            for func_key in self.ovl.case_var_to_fort_var:
                if func_key in d_outputs:
                    func_seeds[func_key] = d_outputs[func_key]
                    if np.abs(func_seeds[func_key]) > 0.0:
                        print(f'  running rev mode derivs for {func_key}')

            csd_seeds = {}
            con_names = self.ovl.get_control_names()
            for func_key in self.ovl.case_derivs_to_fort_var:
                for con_name in con_names:
                    var_name = self.ovl._get_deriv_key(con_name, func_key)

                    if var_name in d_outputs:
                        csd_seeds[var_name] = d_outputs[var_name]
                        
                        if np.abs(csd_seeds[var_name]) > 0.0:
                            # print(var_name, csd_seeds[func_key])
                            print(f'  running rev mode derivs for {var_name}')
                
 
            stab_derivs_seeds = {}
            for func_key in self.ovl.case_stab_derivs_to_fort_var:
                if func_key in d_outputs:
                    stab_derivs_seeds[func_key] = d_outputs[func_key]
                    
                    if np.abs(stab_derivs_seeds[func_key]) > 0.0:
                        # print(var_name, stab_derivs_seeds[func_key])
                        print(f'  running rev mode derivs for {func_key}')

                        
            con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds, gamma_u_seeds, param_seeds, ref_seeds = self.ovl._execute_jac_vec_prod_rev(
                func_seeds=func_seeds, consurf_derivs_seeds=csd_seeds, stab_derivs_seeds=stab_derivs_seeds
            )

            if "gamma" in d_inputs:
                d_inputs["gamma"] += gamma_seeds

            if "gamma_d" in d_inputs:
                d_inputs["gamma_d"] += gamma_d_seeds
            
            if "gamma_u" in d_inputs:
                d_inputs["gamma_u"] += gamma_u_seeds
            
            d_input_geom = om_surf_dict_to_input(geom_seeds)

            for d_input in d_inputs:
                if d_input in d_input_geom:
                    d_inputs[d_input] += d_input_geom[d_input]
                elif d_input in ["alpha", "beta"] or d_input in self.control_names:
                    d_inputs[d_input] += con_seeds[d_input]
                elif d_input in param_seeds:
                    d_inputs[d_input] += param_seeds[d_input]
                elif d_input in ref_seeds:
                    d_inputs[d_input] += ref_seeds[d_input]


# Optional components
class OVLPostProcessComp(om.ExplicitComponent):
    """This component writes out data files for postprocessing. It is optionally added as part of the OVLGroup
    """
    def initialize(self):
        self.options.declare("ovl", types=OVLSolver, recordable=False)
        self.options.declare("output_dir", types=str, recordable=False, default=".")
        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        self.options.declare("write_grid_sol_time", types=bool, default=False)
        
    def setup(self):
        self.ovl = self.options["ovl"]
        self.num_states = self.ovl.get_mesh_size()
        self.num_cs = self.ovl.get_num_control_surfs()
        self.num_vel = self.ovl.NUMAX
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        

        self.add_input("gamma", val=np.zeros(self.num_states))
        self.add_input("gamma_d", val=np.zeros((self.num_cs, self.num_states)))
        self.add_input("gamma_u", val=np.zeros((self.num_vel, self.num_states)))

        add_ovl_conditions_as_inputs(self, self.ovl)
        
        if input_param_vals:
            add_ovl_params_as_inputs(self, self.ovl)
        
        if input_ref_vals:
            add_ovl_refs_as_inputs(self, self.ovl)
            
        add_ovl_controls_as_inputs(self, self.ovl)
        add_ovl_geom_vars(self, self.ovl, add_as="inputs")

        self.res_slice = (slice(0, self.num_states),)
        self.res_d_slice = (slice(0, self.num_cs), slice(0, self.num_states))
        self.res_u_slice = (slice(0, self.num_vel), slice(0, self.num_states))
        
        # check to make sure the output dir exists
        if not os.path.exists(self.options["output_dir"]):
            os.mkdir(self.options["output_dir"])


    def compute(self, inputs, outputs):
        # self.ovl.set_gamma(inputs['gamma'])
        # for c_name in self.control_names:
        #     self.ovl.set_constraint(c_name, inputs[c_name][0])
        # def_dict = self.ovl.get_control_deflections()

        om_set_avl_inputs(self, inputs)

        # update the surface parameters
        surf_data = om_input_to_surf_dict(self, inputs)
        self.ovl.set_surface_params(surf_data)

        gam_arr = inputs["gamma"]
        gam_d_arr = inputs["gamma_d"]
        gam_u_arr = inputs["gamma_u"]

        self.ovl.set_avl_fort_arr("VRTX_R", "GAM", gam_arr, slicer=self.res_slice)

        self.ovl.set_avl_fort_arr("VRTX_R", "GAM_D", gam_d_arr, slicer=self.res_d_slice)
        self.ovl.set_avl_fort_arr("VRTX_R", "GAM_U", gam_u_arr, slicer=self.res_u_slice)

        file_name = f"vlm_{self.iter_count:03d}.avl"
        output_dir = self.options["output_dir"]
        self.ovl.write_geom_file(os.path.join(output_dir, file_name))
        
        file_name = f"vlm_{self.iter_count:03d}"
        if self.options["write_grid_sol_time"]:
            self.ovl.write_tecplot(os.path.join(output_dir, file_name), solution_time=self.iter_count)
        else:
            self.ovl.write_tecplot(os.path.join(output_dir, file_name))
            

class OVLMeshReader(om.ExplicitComponent):
    """
        This class is moslty used to provide an initial set of coordinates for custom paramerization components. It is NOT part of the OVL group
    """
    
    def initialize(self):
        self.options.declare("geom_file", types=str)
        self.options.declare("mass_file", default=None)

    def setup(self):
        geom_file = self.options["geom_file"]
        mass_file = self.options["mass_file"]

        avl = OVLSolver(geo_file=geom_file, mass_file=mass_file, debug=False)
        add_ovl_geom_vars(self, avl, add_as="outputs")


class Differencer(om.ExplicitComponent):
    def setup(self):
        self.add_input('input_vec', shape_by_conn=True)
        # self.add_output('diff_vec',copy_shape='input_vec')
        def compute_shape(shapes):
            # import pdb; pdb.set_trace()
            return (shapes['input_vec'][0]-1,)
        
        self.add_output('diff_vec', compute_shape=compute_shape)
        self.declare_partials("*", "*", method="cs")
        
    def compute(self, inputs, outputs):
        vec = inputs['input_vec']
        diff_vec = vec[1:] - vec[:-1] 
        outputs['diff_vec'] = diff_vec