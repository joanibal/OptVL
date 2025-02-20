# =============================================================================
# Extension modules
# =============================================================================
from optvl import AVLSolver
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
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        # self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        self.avl_solver = AVLSolver(geo_file="aircraft_L1_trans.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.set_constraint("alpha", 5.0)
        self.avl_solver.set_constraint("beta", 0.0)
        self.avl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")


    def finite_dif(self, con_list, geom_seeds, param_seeds, ref_seeds, step=1e-7):
        con_seeds = {}

        for con in con_list:
            con_seeds[con] = 1.0
            
        self.avl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
        self.avl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)
        self.avl_solver.set_parameter_ad_seeds(param_seeds, mode="FD", scale=step)
        self.avl_solver.set_reference_ad_seeds(ref_seeds, mode="FD", scale=step)


        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        # self.avl_solver.execute_run()
        coef_data_peturb = self.avl_solver.get_total_forces()
        consurf_derivs_peturb = self.avl_solver.get_control_stab_derivs()
        stab_deriv_derivs_peturb = self.avl_solver.get_stab_derivs()

        self.avl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-1*step)
        self.avl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-1*step)
        self.avl_solver.set_parameter_ad_seeds(param_seeds, mode="FD", scale=-1*step)
        self.avl_solver.set_reference_ad_seeds(ref_seeds, mode="FD", scale=-1*step)


        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        # self.avl_solver.execute_run()

        coef_data = self.avl_solver.get_total_forces()
        consurf_derivs = self.avl_solver.get_control_stab_derivs()
        stab_deriv_derivs = self.avl_solver.get_stab_derivs()


        func_seeds = {}
        for func_key in coef_data:
            func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

        consurf_derivs_seeds = {}
        for func_key in consurf_derivs:
            consurf_derivs_seeds[func_key] = {}
            for surf_key in consurf_derivs[func_key]:
                consurf_derivs_seeds[func_key][surf_key] = (
                    consurf_derivs_peturb[func_key][surf_key] - consurf_derivs[func_key][surf_key]
                ) / step

        stab_derivs_seeds = {}
        for func_key in stab_deriv_derivs:
            stab_derivs_seeds[func_key] = {}
            for surf_key in stab_deriv_derivs[func_key]:
                stab_derivs_seeds[func_key][surf_key] = (
                    stab_deriv_derivs_peturb[func_key][surf_key] - stab_deriv_derivs[func_key][surf_key]
                ) / step

        return func_seeds, consurf_derivs_seeds, stab_derivs_seeds

    def test_aero_constraint(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.avl_solver.case_var_to_fort_var
        stab_derivs = self.avl_solver.case_stab_derivs_to_fort_var
        sens_funcs = self.avl_solver.execute_run_sensitivies(func_vars)
        sens_sd = self.avl_solver.execute_run_sensitivies([], stab_derivs=stab_derivs, print_timings=False)

        for con_key in self.avl_solver.con_var_to_fort_var:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds = self.finite_dif([con_key], {}, {}, {}, step=1.0e-5)

            for func_key in func_vars:
                ad_dot = sens_funcs[func_key][con_key]
                fd_dot = func_seeds[func_key]

                print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
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

            for func_key in stab_derivs_seeds:
                for var_key in stab_derivs_seeds[func_key]:
                    ad_dot = sens_sd[func_key][var_key][con_key]
                    func_dot = stab_derivs_seeds[func_key][var_key]

                    rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key} wrt {var_key:5}  wrt {con_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    
                    tol = 1e-8
                    if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            atol=1e-9,
                            err_msg=f"{func_key} wrt {var_key:5}  wrt {con_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            rtol=1e-4,
                            err_msg=f"{func_key} wrt {var_key:5}  wrt {con_key}",
                        )
                    
    def test_geom(self):
        # compare the analytical gradients with finite difference for each
        # geometric variable and function

        surf_key = list(self.avl_solver.surf_geom_to_fort_var.keys())[0]
        geom_vars = self.avl_solver.surf_geom_to_fort_var[surf_key]
        cs_names = self.avl_solver.get_control_names()

        consurf_vars = {}
        for func_key in self.avl_solver.case_derivs_to_fort_var:
            consurf_vars[func_key] = [cs_names[0]]
            
        stab_derivs = self.avl_solver.case_stab_derivs_to_fort_var
        # sens = self.avl_solver.execute_run_sensitivies(func_vars)
        func_vars = self.avl_solver.case_var_to_fort_var
        sens = self.avl_solver.execute_run_sensitivies(func_vars, consurf_derivs=consurf_vars, stab_derivs=stab_derivs, print_timings=False)

        # for con_key in self.avl_solver.con_var_to_fort_var:
        sens_FD = {}
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            sens_FD[surf_key] = {}
            for geom_key in geom_vars:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                np.random.seed(arr.size)
                rand_arr = np.random.rand(*arr.shape)
                rand_arr /= np.linalg.norm(rand_arr)

                func_seeds, consurf_deriv_seeds, stab_derivs_seeds = self.finite_dif([], {surf_key: {geom_key: rand_arr}}, {}, {}, step=1.0e-7)

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
                            rtol=1e-4,
                            err_msg=f"{func_key:5} wrt {surf_key}:{geom_key:10}",
                        )

                for func_key in consurf_vars:
                    for cs_key in consurf_vars[func_key]:
                        geom_dot = np.sum(sens[func_key][cs_key][surf_key][geom_key] * rand_arr)
                        func_dot = consurf_deriv_seeds[func_key][cs_key]

                        rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                        # print(
                        #     f"{func_key} wrt {cs_key:5}  wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                        # )
                        
                        tol = 1e-8
                        if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                            # If either value is basically zero, use an absolute tolerance
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                atol=1e-4,
                                err_msg=f"{func_key} wrt {cs_key:5}  wrt {surf_key}:{geom_key:10}",
                            )
                        else:
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                rtol=1e-4,
                                err_msg=f"{func_key} wrt {cs_key:5}  wrt {surf_key}:{geom_key:10}",
                            )
                        
                for func_key in stab_derivs_seeds:
                    for var_key in stab_derivs_seeds[func_key]:
                        geom_dot = np.sum(sens[func_key][var_key][surf_key][geom_key] * rand_arr)
                        func_dot = stab_derivs_seeds[func_key][var_key]


                        rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                        # print(
                        #     f"{func_key} wrt {var_key:5}  wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                        # )
                        
                        tol = 1e-8
                        if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                            # If either value is basically zero, use an absolute tolerance
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                atol=1e-10,
                                err_msg=f"{func_key} wrt {var_key:5}  wrt {surf_key}:{geom_key:10}",
                            )
                        else:
                            np.testing.assert_allclose(
                                geom_dot,
                                func_dot,
                                rtol=1e-3,
                                err_msg=f"{func_key} wrt {var_key:5}  wrt {surf_key}:{geom_key:10}",
                            )                

    def test_params(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.avl_solver.case_var_to_fort_var
        stab_derivs = self.avl_solver.case_stab_derivs_to_fort_var

        sens = self.avl_solver.execute_run_sensitivies(func_vars, stab_derivs=stab_derivs)

        for param_key in self.avl_solver.param_idx_dict:
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
                for var_key in stab_derivs_seeds[func_key]:
                    ad_dot = sens[func_key][var_key][param_key]
                    func_dot = stab_derivs_seeds[func_key][var_key]

                    rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key} wrt {var_key:5}  wrt {param_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    
                    tol = 1e-8
                    if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            atol=1e-10,
                            err_msg=f"{func_key} wrt {var_key:5}  wrt {param_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            rtol=1e-4,
                            err_msg=f"{func_key} wrt {var_key:5}  wrt {param_key}",
                        )

    def test_ref(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.avl_solver.case_var_to_fort_var
        stab_derivs = self.avl_solver.case_stab_derivs_to_fort_var

        sens = self.avl_solver.execute_run_sensitivies(func_vars, stab_derivs=stab_derivs)


        for ref_key in self.avl_solver.ref_var_to_fort_var:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds = self.finite_dif([], {}, {}, {ref_key:1.0}, step=1.0e-5)

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
                for var_key in stab_derivs_seeds[func_key]:
                    ad_dot = sens[func_key][var_key][ref_key]
                    func_dot = stab_derivs_seeds[func_key][var_key]

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
                            err_msg=f"{func_key} wrt {var_key:5}  wrt {ref_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            ad_dot,
                            func_dot,
                            rtol=1e-4,
                            err_msg=f"{func_key} wrt {var_key:5}  wrt {ref_key}",
                        )

if __name__ == "__main__":
    unittest.main()
