# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver

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
geom_dir = os.path.join(base_dir, "..", "geom_files")

geom_file = os.path.join(geom_dir, "aircraft_L1_with_body.avl")
# geom_file = os.path.join(geom_dir, "rect_with_body.avl")
# geom_file = os.path.join(geom_dir, "supra.avl")


class TestTotals(unittest.TestCase):
    # TODO: beta derivatives likely wrong

    def setUp(self):
        self.ovl = OVLSolver(geo_file=geom_file)
        self.ovl.set_variable("alpha", 5.0)
        self.ovl.set_variable("beta", 0.0)
        self.ovl.set_parameter("Mach", 0.8)
        self.ovl.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def finite_dif(self, con_list, geom_seeds, param_seeds, ref_seeds, step=1e-7):
        con_seeds = {}

        for con in con_list:
            con_seeds[con] = 1.0
        self.ovl.set_variable_ad_seeds(con_seeds, mode="FD", scale=step)
        self.ovl.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)
        self.ovl.set_parameter_ad_seeds(param_seeds, mode="FD", scale=step)
        self.ovl.set_reference_ad_seeds(ref_seeds, mode="FD", scale=step)

        self.ovl.avl.update_surfaces()
        self.ovl.avl.get_res()
        self.ovl.avl.exec_rhs()
        self.ovl.avl.get_res()
        self.ovl.avl.velsum()
        self.ovl.avl.aero()
        # self.ovl.execute_run()
        coef_data_peturb = self.ovl.get_total_forces()
        consurf_derivs_peturb = self.ovl.get_control_stab_derivs()
        stab_deriv_derivs_peturb = self.ovl.get_stab_derivs()
        body_axis_deriv_petrub = self.ovl.get_body_axis_derivs()
        body_forces_peturb = self.ovl.get_body_forces()

        self.ovl.set_variable_ad_seeds(con_seeds, mode="FD", scale=-1 * step)
        self.ovl.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-1 * step)
        self.ovl.set_parameter_ad_seeds(param_seeds, mode="FD", scale=-1 * step)
        self.ovl.set_reference_ad_seeds(ref_seeds, mode="FD", scale=-1 * step)

        self.ovl.avl.update_surfaces()
        self.ovl.avl.get_res()
        self.ovl.avl.exec_rhs()
        self.ovl.avl.get_res()
        self.ovl.avl.velsum()
        self.ovl.avl.aero()
        # self.ovl.execute_run()

        coef_data = self.ovl.get_total_forces()
        consurf_derivs = self.ovl.get_control_stab_derivs()
        stab_deriv_derivs = self.ovl.get_stab_derivs()
        body_axis_deriv = self.ovl.get_body_axis_derivs()
        body_forces = self.ovl.get_body_forces()
        
        body_func_seeds = {}
        for body in body_forces:
            body_func_seeds[body] = {}
            for key in body_forces[body]:
                body_func_seeds[body][key] = (body_forces_peturb[body][key] - body_forces[body][key]) / step
        

        func_seeds = {}
        for func_key in coef_data:
            func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

        consurf_derivs_seeds = {}
        for func_key in consurf_derivs:
            consurf_derivs_seeds[func_key] = (consurf_derivs_peturb[func_key] - consurf_derivs[func_key]) / step

        stab_derivs_seeds = {}
        for func_key in stab_deriv_derivs:
            stab_derivs_seeds[func_key] = (stab_deriv_derivs_peturb[func_key] - stab_deriv_derivs[func_key]) / step

        body_axis_derivs_seeds = {}
        for deriv_func in body_axis_deriv:
            body_axis_derivs_seeds[deriv_func] = (
                body_axis_deriv_petrub[deriv_func] - body_axis_deriv[deriv_func]
            ) / step

        return func_seeds, consurf_derivs_seeds, stab_derivs_seeds, body_axis_derivs_seeds

    def test_aero_constraint(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.ovl.case_var_to_fort_var
        stab_derivs = self.ovl.case_stab_derivs_to_fort_var
        body_axis_derivs = self.ovl.case_body_derivs_to_fort_var
        sens_funcs = self.ovl.execute_run_sensitivities(func_vars)
        sens_sd = self.ovl.execute_run_sensitivities([], stab_derivs=stab_derivs, print_timings=False)
        sens_bd = self.ovl.execute_run_sensitivities([], body_axis_derivs=body_axis_derivs, print_timings=False)

        for con_key in self.ovl.con_var_list:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif(
                [con_key], {}, {}, {}, step=1.0e-6
            )

            # for func_key in func_vars:
            for func_key in ['CX']:
                ad_dot = sens_funcs[func_key][con_key]
                fd_dot = func_seeds[func_key]

                # print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))

                # print(f"{func_key:5} wrt {con_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 5e-8
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
                        rtol=5e-4,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

            for func_key in stab_derivs:
                ad_dot = sens_sd[func_key][con_key]
                func_dot = stab_derivs_seeds[func_key]

                rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                # print(
                #     f"{func_key} wrt {con_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                # )

                tol = 5e-8
                if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        atol=5e-8,
                        err_msg=f"{func_key} wrt {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        rtol=5e-4,
                        err_msg=f"{func_key} wrt {con_key}",
                    )

            for func_key in body_axis_derivs_seeds:
                ad_dot = sens_bd[func_key][con_key]
                func_dot = body_axis_derivs_seeds[func_key]

                rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)

                # print(
                #     f"{func_key} wrt {con_key} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                # )

                tol = 5e-8
                if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        atol=2e-8,
                        err_msg=f"{func_key} wrt {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        rtol=1e-3,
                        err_msg=f"{func_key} wrt {con_key}",
                    )

    def test_geom(self):
        # compare the analytical gradients with finite difference for each
        # geometric variable and function

        surf_key = list(self.ovl.surf_geom_to_fort_var.keys())[0]
        geom_vars = self.ovl.surf_geom_to_fort_var[surf_key]
        cs_names = self.ovl.get_control_names()

        consurf_vars = []
        for func_key in self.ovl.case_derivs_to_fort_var:
            consurf_vars.append(self.ovl._get_deriv_key(cs_names[0], func_key))

        func_vars = self.ovl.case_var_to_fort_var
        stab_derivs = self.ovl.case_stab_derivs_to_fort_var
        body_axis_derivs = self.ovl.case_body_derivs_to_fort_var

        sens = self.ovl.execute_run_sensitivities(
            func_vars,
            consurf_derivs=consurf_vars,
            stab_derivs=stab_derivs,
            body_axis_derivs=body_axis_derivs,
            print_timings=False,
        )

        # for con_key in self.ovl.con_var_list:
        sens_FD = {}
        for surf_key in self.ovl.surf_geom_to_fort_var:
            sens_FD[surf_key] = {}
            for geom_key in geom_vars:
                arr = self.ovl.get_surface_param(surf_key, geom_key)
                np.random.seed(arr.size)
                rand_arr = np.random.rand(*arr.shape)
                rand_arr /= np.linalg.norm(rand_arr)

                func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif(
                    [], {surf_key: {geom_key: rand_arr}}, {}, {}, step=1.0e-7
                )

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

                    # rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)
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
                            rtol=6e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )

                for func_key in stab_derivs_seeds:
                    geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                    func_dot = stab_derivs_seeds[func_key]

                    rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key}  wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )

                    tol = 5e-7
                    if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            atol=5e-9,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            rtol=6e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )

                for func_key in body_axis_derivs_seeds:
                    geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                    func_dot = body_axis_derivs_seeds[func_key]

                    rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key}  wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )

                    tol = 1e-6
                    if np.abs(geom_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            atol=5e-8,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            geom_dot,
                            func_dot,
                            rtol=6e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{geom_key:10}",
                        )

    def test_params(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.ovl.case_var_to_fort_var
        stab_derivs = self.ovl.case_stab_derivs_to_fort_var

        sens = self.ovl.execute_run_sensitivities(func_vars, stab_derivs=stab_derivs)

        for param_key in self.ovl.param_idx_dict:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif(
                [], {}, {param_key: 1.0}, {}, step=1.0e-6
            )

            for func_key in func_vars:
                ad_dot = sens[func_key][param_key]
                fd_dot = func_seeds[func_key]

                # rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))
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

                # rel_err = np.abs(ad_dot - func_dot) / np.abs(func_dot + 1e-20)
                # print(
                #     f"{func_key:20} wrt {param_key:10} | AD:{ad_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                # )

                tol = 1e-8
                if np.abs(ad_dot) < tol or np.abs(func_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        atol=1e-9,
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
        func_vars = self.ovl.case_var_to_fort_var
        stab_derivs = self.ovl.case_stab_derivs_to_fort_var

        sens = self.ovl.execute_run_sensitivities(func_vars, stab_derivs=stab_derivs)

        for ref_key in self.ovl.ref_var_to_fort_var:
            # for con_key in ['beta']:
            func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif(
                [], {}, {}, {ref_key: 1.0}, step=1.0e-5
            )

            for func_key in func_vars:
                ad_dot = sens[func_key][ref_key]
                fd_dot = func_seeds[func_key]

                # print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))

                # print(f"{func_key:5} wrt {ref_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-13
                if np.abs(np.linalg.norm(ad_dot)) < tol or np.abs(fd_dot) < tol:
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
                if np.abs(np.linalg.norm(ad_dot)) < tol or np.abs(func_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        atol=1e-9,
                        err_msg=f"{func_key}  wrt {ref_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        func_dot,
                        rtol=1e-4,
                        err_msg=f"{func_key}  wrt {ref_key}",
                    )

class TestDirectVsAdjoint(unittest.TestCase):

    def setUp(self):
        self.ovl = OVLSolver(geo_file=geom_file)
        self.ovl.set_variable("alpha", 5.0)
        self.ovl.set_variable("beta", 0.0)
        self.ovl.set_parameter("Mach", 0.8)
        self.ovl.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_aero_constraint(self):
        # compare the analytical gradients with finite difference for each constraint and function
        # func_vars = self.ovl.case_var_to_fort_var
        stab_derivs = ["dCL/dalpha"]
        # body_axis_derivs = self.ovl.case_body_derivs_to_fort_var
        funcs = ["CL", "CD", "Cm"]
        con_dvs = ["alpha"]
        ref_dvs = ["Sref"]
        param_dvs = ["Mach"]
        geom_dvs = [("Wing", "scale"),("Wing", "chords")]
        
        sens_adjoint = self.ovl.execute_run_sensitivities(funcs, stab_derivs=stab_derivs )
        
        sens_direct = self.ovl.execute_run_sensitivities_direct(geom_dvs=geom_dvs, con_dvs=con_dvs, ref_dvs=ref_dvs, param_dvs=param_dvs,  add_stab_derivs=True)
        
        for func in funcs + stab_derivs:
            for dv in geom_dvs:
                # print(f"Adjoint d{func}/d{[dv[0]]} {[dv[1]]} {sens_adjoint[func][dv[0]][dv[1]]}")
                # print(f"Direct  d{func}/d{[dv[0]]} {[dv[1]]} {sens_direct[func][dv[0]][dv[1]]}")
                np.testing.assert_allclose(sens_direct[func][dv[0]][dv[1]], sens_adjoint[func][dv[0]][dv[1]], 1e-15, 1e-15)
                
            
            for dv in con_dvs + ref_dvs + param_dvs:
                # print(f"Adjoint d{func}/d{dv} {sens_adjoint[func][dv]}")
                # print(f"Direct  d{func}/d{dv} {sens_direct[func][dv]}")
                np.testing.assert_allclose(sens_direct[func][dv], sens_adjoint[func][dv], 1e-15, 1e-15)
            

if __name__ == "__main__":
    unittest.main()
