# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver
import copy
import platform

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


class TestFunctionPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        # self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl")
        self.ovl_solver = OVLSolver(geo_file="rect.avl")
        self.ovl_solver.set_constraint("alpha", 25.0)
        self.ovl_solver.set_constraint("beta", 5.0)
        self.ovl_solver.execute_run()
        
    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")


    def test_fwd_aero_constraint(self):
        for con_key in self.ovl_solver.con_var_to_fort_var:
            func_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})[0]

            func_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, geom_seeds={}, mode="FD", step=1e-7
            )[0]

            for func_key in func_seeds:
                # print(f"{func_key} wrt {con_key}", func_seeds[func_key], func_seeds_FD[func_key])
                tol = 1e-13
                if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        atol=1e-6,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        rtol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

    def test_rev_aero_constraint(self):
        for con_key in self.ovl_solver.con_var_to_fort_var:
            self.ovl_solver.clear_ad_seeds_fast()

            func_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})[0]
            self.ovl_solver.clear_ad_seeds_fast()

            for func_key in self.ovl_solver.case_var_to_fort_var:
                con_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})[0]
                # print(f"{func_key} wrt {con_key}", "fwd ", func_seeds_fwd[func_key], "rev", con_seeds_rev[con_key])
                tol = 1e-14

                if np.abs(func_seeds_fwd[func_key]) < tol or np.abs(con_seeds_rev[con_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        con_seeds_rev[con_key],
                        atol=1e-14,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        con_seeds_rev[con_key],
                        rtol=1e-12,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                func_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[0]

                func_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-7
                )[0]

                for func_key in func_seeds:
                    rel_error = np.linalg.norm(func_seeds[func_key] - func_seeds_FD[func_key]) / np.linalg.norm(
                        func_seeds_FD[func_key] + 1e-15
                    )

                    # print(
                    #     f"{func_key:10} wrt {surf_key}:{geom_key} AD:{func_seeds[func_key]: .5e} FD:{func_seeds_FD[func_key]: .5e} rel_error:{rel_error: .3e}"
                    # )

                    tol = 1e-13
                    if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            func_seeds[func_key],
                            func_seeds_FD[func_key],
                            atol=1e-6,
                            err_msg=f"func_key {func_key} w.r.t. {geom_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            func_seeds[func_key],
                            func_seeds_FD[func_key],
                            rtol=5e-4,
                            err_msg=f"func_key {func_key} w.r.t. {geom_key}",
                        )

    def test_rev_geom(self):
        np.random.seed(111)

        sens_dict_rev = {}
        for func_key in self.ovl_solver.case_var_to_fort_var:
            sens_dict_rev[func_key] = self.ovl_solver._execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})[1]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                func_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[0]

                for func_key in func_seeds_fwd:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(sens_dict_rev[func_key][surf_key][geom_key] * geom_seeds)
                    fwd_sum = np.sum(func_seeds_fwd[func_key])

                    # # print(geom_seeds_rev)
                    tol = 1e-13
                    # print(f"{func_key} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(fwd_sum) < tol or np.abs(rev_sum) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            atol=1e-14,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{geom_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            rtol=1e-12,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{geom_key}",
                        )

    def test_fwd_gamma(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        np.random.seed(2)
        gamma_seeds = np.random.rand(num_gamma)

        func_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds)[0]
        func_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds, mode="FD", step=1e-7)[0]

        for func_key in func_seeds:
            rel_error = np.linalg.norm(func_seeds[func_key] - func_seeds_FD[func_key]) / np.linalg.norm(
                func_seeds_FD[func_key] + 1e-15
            )
            # print(
            #     f"{func_key:10} AD:{func_seeds[func_key]: .5e} FD:{func_seeds_FD[func_key]: .5e} rel_error:{rel_error: .3e}"
            # )

            tol = 1e-13
            if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                # If either value is basically zero, use an absolute tolerance
                np.testing.assert_allclose(
                    func_seeds[func_key],
                    func_seeds_FD[func_key],
                    atol=1e-6,
                    err_msg=f"func_key {func_key} w.r.t. gamma",
                )
            else:
                np.testing.assert_allclose(
                    func_seeds[func_key],
                    func_seeds_FD[func_key],
                    rtol=5e-5,
                    err_msg=f"func_key {func_key} w.r.t. gamma",
                )

    def test_rev_gamma(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        gamma_seeds_fwd = np.random.rand(num_gamma)

        func_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds_fwd)[0]
        self.ovl_solver.clear_ad_seeds_fast()

        for func_key in self.ovl_solver.case_var_to_fort_var:
            gamma_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})[2]

            rev_sum = np.sum(gamma_seeds_rev * gamma_seeds_fwd)
            fwd_sum = np.sum(func_seeds_fwd[func_key])

            np.testing.assert_allclose(
                fwd_sum,
                rev_sum,
                atol=1e-14,
                err_msg=f"func_key {func_key} w.r.t. gamma",
            )


    def test_fwd_param(self):
        for param_key in self.ovl_solver.param_idx_dict:
            func_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(param_seeds={param_key: 1.0})[0]

            func_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                param_seeds={param_key: 1.0}, mode="FD", step=1e-7
            )[0]

            for func_key in func_seeds:
                # print(f"{func_key} wrt {param_key}", func_seeds[func_key], func_seeds_FD[func_key])
                tol = 1e-13
                if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        atol=1e-6,
                        err_msg=f"func_key {func_key} w.r.t. {param_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        rtol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {param_key}",
                    )

    def test_rev_param(self):
        for param_key in self.ovl_solver.param_idx_dict:
            self.ovl_solver.clear_ad_seeds_fast()

            func_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(param_seeds={param_key: 1.0})[0]
            self.ovl_solver.clear_ad_seeds_fast()

            for func_key in self.ovl_solver.case_var_to_fort_var:
                param_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})[5]
                # print(f"{func_key} wrt {param_key}", "fwd ", func_seeds_fwd[func_key], "rev", param_seeds_rev[param_key])
                tol = 1e-14

                if np.abs(func_seeds_fwd[func_key]) < tol or np.abs(param_seeds_rev[param_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        param_seeds_rev[param_key],
                        atol=1e-14,
                        err_msg=f"func_key {func_key} w.r.t. {param_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        param_seeds_rev[param_key],
                        rtol=1e-12,
                        err_msg=f"func_key {func_key} w.r.t. {param_key}",
                    )


    def test_fwd_ref(self):
        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            func_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[0]

            func_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                ref_seeds={ref_key: 1.0}, mode="FD", step=1e-7
            )[0]

            for func_key in func_seeds:
                # print(f"{func_key} wrt {ref_key}", func_seeds[func_key], func_seeds_FD[func_key])
                tol = 1e-13
                if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        atol=1e-6,
                        err_msg=f"func_key {func_key} w.r.t. {ref_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        rtol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {ref_key}",
                    )

    def test_rev_ref(self):
        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            self.ovl_solver.clear_ad_seeds_fast()

            func_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[0]
            self.ovl_solver.clear_ad_seeds_fast()

            for func_key in self.ovl_solver.case_var_to_fort_var:
                ref_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})[6]
                
                # print(f"{func_key} wrt {ref_key}", "fwd ", func_seeds_fwd[func_key], "rev", ref_seeds_rev[ref_key])
                tol = 1e-14

                if np.abs(func_seeds_fwd[func_key]) < tol or np.abs(ref_seeds_rev[ref_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        ref_seeds_rev[ref_key],
                        atol=1e-14,
                        err_msg=f"func_key {func_key} w.r.t. {ref_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        ref_seeds_rev[ref_key],
                        rtol=1e-12,
                        err_msg=f"func_key {func_key} w.r.t. {ref_key}",
                    )


class TestResidualPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl")
        # self.ovl_solver = OVLSolver(geo_file="rect.avl")
        self.ovl_solver.set_constraint("alpha", 25.0)
        self.ovl_solver.set_constraint("beta", 5.0)
        self.ovl_solver.execute_run()
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    
    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")


    def test_fwd_aero_constraint(self):
        for con_key in self.ovl_solver.con_var_to_fort_var:
            res_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, geom_seeds={}, mode="FD", step=1e-8
            )[1]
            
            res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})[1]



            # print(f"res wrt {con_key}", np.linalg.norm(res_seeds), np.linalg.norm(res_seeds_FD))
            np.testing.assert_allclose(
                res_seeds,
                res_seeds_FD,
                rtol=1e-5,
                err_msg=f"d(res) w.r.t.{con_key}"
            )

    def test_rev_aero_constraint(self):
        num_res = self.ovl_solver.get_mesh_size()
        res_seeds_rev = np.random.rand(num_res)
        con_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)[0]

        self.ovl_solver.clear_ad_seeds_fast()

        for con_key in self.ovl_solver.con_var_to_fort_var:
            res_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[1]

            # do dot product
            res_sum = np.sum(res_seeds_rev * res_seeds_fwd)
            con_sum = np.sum(con_seeds_rev[con_key])

            # print(f"res wrt {con_key}", "rev", con_sum, "fwd", res_sum)

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[1]

                res_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )[1]

                abs_error = np.abs(res_seeds - res_seeds_FD)
                rel_error = (res_seeds - res_seeds_FD) / (res_seeds + 1e-15)
                idx_max_rel_error = np.argmax(np.abs(rel_error))
                idx_max_abs_error = np.argmax(np.abs(abs_error))

                # print(
                #     f"{surf_key:10} {geom_key:10} AD:{np.linalg.norm(res_seeds): .5e} FD:{np.linalg.norm(res_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                # )
                np.testing.assert_allclose(
                    res_seeds,
                    res_seeds_FD,
                    atol=1e-5,
                    err_msg=f"func_key res w.r.t. {surf_key}:{geom_key}",
                )

    def test_rev_geom(self):
        num_res = self.ovl_solver.get_mesh_size()
        res_seeds_rev = np.random.seed(111)
        res_seeds_rev = np.random.rand(num_res)
        geom_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)[1]

        self.ovl_solver.clear_ad_seeds_fast()
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[1]

                # do dot product
                res_sum = np.sum(res_seeds_rev * res_seeds)
                geom_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds)

                # print(f"res wrt {surf_key}:{geom_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    geom_sum,
                    atol=1e-14,
                    err_msg=f"func_key res w.r.t. {surf_key}:{geom_key}",
                )

    def test_fwd_gamma(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        np.random.seed(111)
        gamma_seeds = np.random.rand(num_gamma)

        res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds)[1]
        res_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds, mode="FD", step=1e-7)[1]

        np.testing.assert_allclose(
            res_seeds,
            res_seeds_FD,
            rtol=1e-5,
        )


    def test_fwd_param(self):
        for param_key in self.ovl_solver.param_idx_dict:
            res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(param_seeds={param_key: 1.0})[1]

            res_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                param_seeds={param_key: 1.0}, mode="FD", step=1e-7
            )[1]

            # print(f"res wrt {param_key}", np.linalg.norm(res_seeds), np.linalg.norm(res_seeds_FD))
            np.testing.assert_allclose(
                res_seeds,
                res_seeds_FD,
                atol=1e-6,
                err_msg=f"d(res) w.r.t.{param_key}"
            )

    def test_rev_param(self):
        num_res = self.ovl_solver.get_mesh_size()
        res_seeds_rev = np.random.rand(num_res)
        param_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)[5]

        self.ovl_solver.clear_ad_seeds_fast()

        for param_key in self.ovl_solver.param_idx_dict:
            res_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(param_seeds={param_key: 1.0})[1]

            # do dot product
            res_sum = np.sum(res_seeds_rev * res_seeds_fwd)
            param_sum = np.sum(param_seeds_rev[param_key])

            # print(f"res wrt {param_key}", "rev", param_sum, "fwd", res_sum)
            np.testing.assert_allclose(
                res_sum,
                param_sum,
                atol=1e-14,
                err_msg=f"func_key res w.r.t. {param_key}",
            )
            

    def test_fwd_ref(self):
        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[1]
            res_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                ref_seeds={ref_key: 1.0}, mode="FD", step=1e-7
            )[1]

            # print(f"res wrt {ref_key}", np.linalg.norm(res_seeds), np.linalg.norm(res_seeds_FD))
            np.testing.assert_allclose(
                res_seeds,
                res_seeds_FD,
                atol=1e-6,
                err_msg=f"d(res) w.r.t.{ref_key}"
            )

    def test_rev_ref(self):
        num_res = self.ovl_solver.get_mesh_size()
        res_seeds_rev = np.random.rand(num_res)
        ref_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)[6]

        self.ovl_solver.clear_ad_seeds_fast()

        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            res_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[1]
            # do dot product
            res_sum = np.sum(res_seeds_rev * res_seeds_fwd)
            ref_sum = np.sum(ref_seeds_rev[ref_key])

            # print(f"res wrt {ref_key}", "rev", ref_sum, "fwd", res_sum)
            np.testing.assert_allclose(
                res_sum,
                ref_sum,
                atol=1e-14,
                err_msg=f"func_key res w.r.t. {ref_key}",
            )
            




if __name__ == "__main__":
    unittest.main()
