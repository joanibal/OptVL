# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver
import copy

# =============================================================================
# Standard Python Modules
# =============================================================================
import os
import psutil
# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestTotals(unittest.TestCase):
    # TODO: beta derivatives likely wrong

    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        # self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl")
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1_trans.avl")
        # self.ovl_solver = OVLSolver(geo_file="rect.avl")
        self.ovl_solver.set_constraint("alpha", 5.0)
        self.ovl_solver.set_constraint("beta", 0.0)
        self.ovl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")


    def finite_dif(self, con_list, geom_seeds, param_seeds, ref_seeds, step=1e-7):
        con_seeds = {}

        for con in con_list:
            con_seeds[con] = 1.0
            
        self.ovl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
        self.ovl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)
        self.ovl_solver.set_parameter_ad_seeds(param_seeds, mode="FD", scale=step)
        self.ovl_solver.set_reference_ad_seeds(ref_seeds, mode="FD", scale=step)


        self.ovl_solver.avl.update_surfaces()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.exec_rhs()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.velsum()
        self.ovl_solver.avl.aero()
        # self.ovl_solver.execute_run()
        coef_data_peturb = self.ovl_solver.get_total_forces()
        consurf_derivs_peturb = self.ovl_solver.get_control_stab_derivs()
        stab_deriv_derivs_peturb = self.ovl_solver.get_stab_derivs()
        body_axis_deriv_petrub = self.ovl_solver.get_body_axis_derivs()

        self.ovl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-1*step)
        self.ovl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-1*step)
        self.ovl_solver.set_parameter_ad_seeds(param_seeds, mode="FD", scale=-1*step)
        self.ovl_solver.set_reference_ad_seeds(ref_seeds, mode="FD", scale=-1*step)


        self.ovl_solver.avl.update_surfaces()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.exec_rhs()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.velsum()
        self.ovl_solver.avl.aero()
        # self.ovl_solver.execute_run()

        coef_data = self.ovl_solver.get_total_forces()
        consurf_derivs = self.ovl_solver.get_control_stab_derivs()
        stab_deriv_derivs = self.ovl_solver.get_stab_derivs()
        body_axis_deriv = self.ovl_solver.get_body_axis_derivs()

        func_seeds = {}
        for func_key in coef_data:
            func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

        consurf_derivs_seeds = {}
        for func_key in consurf_derivs:
            consurf_derivs_seeds[func_key] = (
                consurf_derivs_peturb[func_key] - consurf_derivs[func_key]
            ) / step

        stab_derivs_seeds = {}
        for func_key in stab_deriv_derivs:
            stab_derivs_seeds[func_key] = (
                stab_deriv_derivs_peturb[func_key] - stab_deriv_derivs[func_key]
            ) / step

        body_axis_derivs_seeds = {}
        for deriv_func in body_axis_deriv:
            body_axis_derivs_seeds[deriv_func] = (
                body_axis_deriv_petrub[deriv_func] - body_axis_deriv[deriv_func]
            ) / step

        return func_seeds, consurf_derivs_seeds, stab_derivs_seeds, body_axis_derivs_seeds

    def test_aero_constraint(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.ovl_solver.case_var_to_fort_var
        stab_derivs = self.ovl_solver.case_stab_derivs_to_fort_var
        body_axis_derivs = self.ovl_solver.case_body_derivs_to_fort_var
        sens_funcs = self.ovl_solver.execute_run_sensitivities(func_vars)
        sens_sd = self.ovl_solver.execute_run_sensitivities([], stab_derivs=stab_derivs, print_timings=False)
        sens_bd = self.ovl_solver.execute_run_sensitivities([], body_axis_derivs=body_axis_derivs, print_timings=False)

        for con_key in self.ovl_solver.con_var_to_fort_var:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif([con_key], {}, {}, {}, step=1.0e-5)

            for func_key in func_vars:
                ad_dot = sens_funcs[func_key][con_key]
                fd_dot = func_seeds[func_key]

                # print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))

                # print(f"{func_key:5} wrt {con_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-8
                if np.abs(ad_dot) < tol or np.abs(fd_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        atol=1e-9,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        rtol=5e-5,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

            for func_key in stab_derivs:
                    ad_dot = sens_sd[func_key][con_key]
                    func_dot = stab_derivs_seeds[func_key]

                    rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key} wrt {con_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    
                    tol = 1e-8
                    if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            atol=1e-9,
                            err_msg=f"{func_key} wrt {con_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            rtol=1e-4,
                            err_msg=f"{func_key} wrt {con_key}",
                        )
            
            for func_key in body_axis_derivs_seeds:
                    ad_dot = sens_bd[func_key][con_key]
                    func_dot = body_axis_derivs_seeds[func_key]

                    rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                    print(
                        f"{func_key} wrt {con_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    )
                    
                    tol = 1e-8
                    if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            atol=1e-9,
                            err_msg=f"{func_key} wrt {con_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            rtol=1e-4,
                            err_msg=f"{func_key} wrt {con_key}",
                        )
                    
    def test_geom(self):
        # compare the analytical gradients with finite difference for each
        # geometric variable and function

        surf_key = list(self.ovl_solver.surf_geom_to_fort_var.keys())[0]
        geom_vars = self.ovl_solver.surf_geom_to_fort_var[surf_key]
        cs_names = self.ovl_solver.get_control_names()

        consurf_vars = []
        for func_key in self.ovl_solver.case_derivs_to_fort_var:
            consurf_vars.append(self.ovl_solver._get_deriv_key(cs_names[0], func_key))
   
            
        func_vars = self.ovl_solver.case_var_to_fort_var
        stab_derivs = self.ovl_solver.case_stab_derivs_to_fort_var
        body_axis_derivs = self.ovl_solver.case_body_derivs_to_fort_var
        
        sens = self.ovl_solver.execute_run_sensitivities(func_vars, consurf_derivs=consurf_vars, stab_derivs=stab_derivs, body_axis_derivs=body_axis_derivs, print_timings=False)

        # for con_key in self.ovl_solver.con_var_to_fort_var:
        sens_FD = {}
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            sens_FD[surf_key] = {}
            for geom_key in geom_vars:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                np.random.seed(arr.size)
                rand_arr = np.random.rand(*arr.shape)
                rand_arr /= np.linalg.norm(rand_arr)

                func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif([], {surf_key: {geom_key: rand_arr}}, {}, {}, step=1.0e-7)

                for func_key in func_vars:
                    geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                    func_dot = func_seeds[func_key]

                    rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key:5} wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    tol = 1e-7
                    if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            atol=1e-4,
                            err_msg=f"{func_key:5} wrt {surf_key}:{geom_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            rtol=5e-3,
                            err_msg=f"{func_key:5} wrt {surf_key}:{geom_key:10}",
                        )

                for func_key in consurf_vars:
                    # for cs_key in consurf_vars[func_key]:
                    geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                    func_dot = consurf_deriv_seeds[func_key]

                    rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key} wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    
                    tol = 1e-8
                    if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            atol=1e-4,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            rtol=3e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )
                        
                for func_key in stab_derivs_seeds:
                        geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                        func_dot = stab_derivs_seeds[func_key]


                        rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                        # print(
                        #     f"{func_key}  wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                        # )
                        
                        tol = 1e-10
                        if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                            # If either value is basically zero, use an absolute tolerance
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                atol=1e-7,
                                err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                            )
                        else:
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                rtol=3e-3,
                                err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                            )                
                
                for func_key in body_axis_derivs_seeds:
                        geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                        func_dot = body_axis_derivs_seeds[func_key]


                        rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                        print(
                            f"{func_key}  wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                        )
                        
                        tol = 5e-7
                        if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                            # If either value is basically zero, use an absolute tolerance
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                atol=tol,
                                err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                            )
                        else:
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                rtol=3e-3,
                                err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                            )                
                
                    

    def test_params(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.ovl_solver.case_var_to_fort_var
        stab_derivs = self.ovl_solver.case_stab_derivs_to_fort_var

        sens = self.ovl_solver.execute_run_sensitivities(func_vars, stab_derivs=stab_derivs)

        for param_key in self.ovl_solver.param_idx_dict:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds  = self.finite_dif([], {}, {param_key:1.0}, {}, step=1.0e-6)

            for func_key in func_vars:
                ad_dot = sens[func_key][param_key]
                fd_dot = func_seeds[func_key]

                # print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))

                # print(f"{func_key:5} wrt {param_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-13
                if np.abs(ad_dot) < tol or np.abs(fd_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        atol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {param_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        rtol=5e-4,
                        err_msg=f"func_key {func_key} w.r.t. {param_key}",
                    )

            for func_key in stab_derivs_seeds:
                    ad_dot = sens[func_key][param_key]
                    func_dot = stab_derivs_seeds[func_key]

                    rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key}   wrt {param_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    
                    tol = 1e-8
                    if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            atol=1e-10,
                            err_msg=f"{func_key}  wrt {param_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            rtol=1e-4,
                            err_msg=f"{func_key}  wrt {param_key}",
                        )

    def test_ref(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.ovl_solver.case_var_to_fort_var
        stab_derivs = self.ovl_solver.case_stab_derivs_to_fort_var

        sens = self.ovl_solver.execute_run_sensitivities(func_vars, stab_derivs=stab_derivs)


        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif([], {}, {}, {ref_key:1.0}, step=1.0e-5)

            for func_key in func_vars:
                ad_dot = sens[func_key][ref_key]
                fd_dot = func_seeds[func_key]

                # print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))

                # print(f"{func_key:5} wrt {ref_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-13
                if np.abs(ad_dot) < tol or np.abs(fd_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        atol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {ref_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        rtol=5e-4,
                        err_msg=f"func_key {func_key} w.r.t. {ref_key}",
                    )
            
            for func_key in stab_derivs_seeds:
                ad_dot = sens[func_key][ref_key]
                func_dot = stab_derivs_seeds[func_key]

                # rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                # print(
                #     f"{func_key} wrt {var_key:5}  wrt {ref_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                # )
                tol = 1e-8
                if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        atol=1e-10,
                        err_msg=f"{func_key}  wrt {ref_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        rtol=1e-4,
                        err_msg=f"{func_key}  wrt {ref_key}",
                    )

if __name__ == "__main__":
    unittest.main()
