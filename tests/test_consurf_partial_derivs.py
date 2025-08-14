# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver
import copy
import psutil

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestResidualDPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl")
        # self.ovl_solver = OVLSolver(geo_file="rect.avl")
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
            res_d_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[5]

            res_d_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, mode="FD", step=1e-5
            )[5]

            np.testing.assert_allclose(
                res_d_seeds,
                res_d_seeds_FD,
                rtol=1e-5,
            )

    def test_rev_aero_constraint(self):
        for con_key in self.ovl_solver.con_var_to_fort_var:
            res_d_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[5]

            num_gamma = self.ovl_solver.get_mesh_size()
            num_consurf = self.ovl_solver.get_num_control_surfs()
            res_d_seeds_rev = np.random.rand(num_consurf, num_gamma)
            res_d_seeds_rev = np.ones_like(res_d_seeds_rev)

            self.ovl_solver.clear_ad_seeds_fast()

            con_seeds = self.ovl_solver._execute_jac_vec_prod_rev(res_d_seeds=res_d_seeds_rev)[0]
            self.ovl_solver.clear_ad_seeds_fast()

            # do dot product
            res_sum = np.sum(res_d_seeds_rev * res_d_seeds)

            # print(f"res wrt {con_key}", "fwd", res_sum, "rev", con_seeds[con_key])

            np.testing.assert_allclose(
                res_sum,
                con_seeds[con_key],
                atol=1e-14,
                err_msg=f"deriv_func res w.r.t. {con_key}",
            )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                res_d_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[5]

                res_d_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )[5]
                abs_error = np.abs(res_d_seeds.flatten() - res_d_seeds_FD.flatten())
                rel_error = np.abs((res_d_seeds.flatten() - res_d_seeds_FD.flatten()) / (res_d_seeds.flatten() + 1e-15))

                idx_max_rel_error = np.argmax(rel_error)
                idx_max_abs_error = np.argmax(abs_error)
                # print(
                #     f"{surf_key:10} {geom_key:10} AD:{np.linalg.norm(res_d_seeds): .5e} FD:{np.linalg.norm(res_d_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                # )
                np.testing.assert_allclose(
                    res_d_seeds,
                    res_d_seeds_FD,
                    atol=1e-4,
                    err_msg=f"deriv_func res w.r.t. {surf_key}:{geom_key}",
                )

    def test_rev_geom(self):
        np.random.seed(111)
        num_gamma = self.ovl_solver.get_mesh_size()
        num_consurf = self.ovl_solver.get_num_control_surfs()
        res_d_seeds_rev = np.random.rand(num_consurf, num_gamma)
        res_d_seeds_rev = np.ones_like(res_d_seeds_rev)

        self.ovl_solver.clear_ad_seeds_fast()

        geom_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_d_seeds=res_d_seeds_rev)[1]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                res_d_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[5]

                # do dot product
                res_sum = np.sum(res_d_seeds_rev * res_d_seeds_fwd)
                geom_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds)

                # print(f"res wrt {surf_key}:{geom_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    geom_sum,
                    atol=1e-14,
                    err_msg=f"deriv_func res w.r.t. {surf_key}:{geom_key}",
                )
        self.ovl_solver.clear_ad_seeds_fast()

    def test_fwd_gamma_d(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        num_consurf = self.ovl_solver.get_num_control_surfs()
        gamma_d_seeds = np.random.rand(num_consurf, num_gamma)
        # gamma_d_seeds = np.array([[1],[0]])

        res_d_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds)[5]

        res_d_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
            gamma_d_seeds=gamma_d_seeds, mode="FD", step=1e-0
        )[5]

        np.testing.assert_allclose(
            res_d_seeds,
            res_d_seeds_FD,
            rtol=1e-5,
        )

    def test_rev_gamma_d(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        num_consurf = self.ovl_solver.get_num_control_surfs()
        gamma_d_seeds_fwd = np.random.rand(num_consurf, num_gamma)

        res_d_seeds_rev = np.random.rand(num_consurf, num_gamma)

        res_d_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds_fwd)[5]
        self.ovl_solver.clear_ad_seeds_fast()

        gamma_d_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_d_seeds=res_d_seeds_rev)[3]

        gamma_sum = np.sum(gamma_d_seeds_rev * gamma_d_seeds_fwd)
        res_sum = np.sum(res_d_seeds_rev * res_d_seeds_fwd)

        # print("fwd_sum", gamma_sum, "rev_sum", res_sum)
        np.testing.assert_allclose(
            gamma_sum,
            res_sum,
            atol=1e-14,
            err_msg=f"res w.r.t. gamma",
        )


class TestConSurfDerivsPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl")
        # self.ovl_solver = OVLSolver(geo_file="rect.avl")
        self.ovl_solver.set_constraint("alpha", 45.0)
        self.ovl_solver.set_constraint("beta", 45.0)
        self.ovl_solver.execute_run()
        self.ovl_solver.clear_ad_seeds_fast()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_aero_constraint(self):
        for con_key in self.ovl_solver.con_var_to_fort_var:
            cs_d = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[2]

            cs_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, mode="FD", step=1e-8)[2]
            
            for deriv_func in cs_d:
                sens_label = f"d{deriv_func} wrt {con_key}"
                np.testing.assert_allclose(
                    cs_d[deriv_func],
                    cs_d_fd[deriv_func],
                    rtol=1e-3,
                    err_msg=sens_label,
                )

    def test_rev_aero_constraint(self):
        cs_names = self.ovl_solver.get_control_names()

        cs_deriv_seeds = {}
        for deriv_func in self.ovl_solver.case_derivs_to_fort_var:
            for cs_name in cs_names:
                cs_deriv_seeds[f'd{deriv_func}/d{cs_name}'] = np.random.rand(1)[0]

        con_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_deriv_seeds)[0]

        self.ovl_solver.clear_ad_seeds_fast()

        for con_key in self.ovl_solver.con_var_to_fort_var:
            cs_deriv_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[2]
 
            cs_deriv_sum = 0.0
            for deriv_func in cs_deriv_seeds_fwd:
                cs_deriv_sum += cs_deriv_seeds[deriv_func] * cs_deriv_seeds_fwd[deriv_func]

            # do dot product
            con_sum = np.sum(con_seeds_rev[con_key])

            # print(f"cs_dervs wrt {con_key}", "rev", con_sum, "fwd", cs_deriv_sum)

            np.testing.assert_allclose(
                con_sum,
                cs_deriv_sum,
                atol=1e-14,
                err_msg=f"cs_dervs wrt {con_key}",
            )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                cs_d = self.ovl_solver._execute_jac_vec_prod_fwd(geom_seeds={surf_key: {geom_key: geom_seeds}})[2]

                cs_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )[2]

                for deriv_func in cs_d:
                    sens_label = f"{deriv_func} wrt {surf_key}:{geom_key:5}"

                    # print(f"{sens_label} AD:{cs_d[deriv_func]} FD:{cs_d_fd[deriv_func]}")

                    tol = 1e-10
                    # print(f"{deriv_func} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(cs_d[deriv_func]) < tol or np.abs(cs_d_fd[deriv_func]) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            cs_d[deriv_func],
                            cs_d_fd[deriv_func],
                            atol=5e-13,
                            err_msg=sens_label,
                        )
                    else:
                        np.testing.assert_allclose(
                            cs_d[deriv_func],
                            cs_d_fd[deriv_func],
                            rtol=1e-3,
                            err_msg=sens_label,
                            )

    def test_rev_geom(self):
        np.random.seed(111)
        cs_names = self.ovl_solver.get_control_names()

        cs_d_rev = {}
        cs_d_rev = {}
        for deriv_func in self.ovl_solver.case_derivs_to_fort_var:
            for cs_name in cs_names:
                cs_d_rev[f'd{deriv_func}/d{cs_name}'] = np.random.rand(1)[0]

        geom_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_d_rev)[1]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds_fwd = np.random.rand(*arr.shape)

                func_seeds_fwd, _, cs_d_fwd, _, _, _,_ = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds_fwd}}
                )
                

                for deriv_func in func_seeds_fwd:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds_fwd)

                    fwd_sum = 0.0
                    for deriv_func in cs_d_fwd:
                        fwd_sum += cs_d_rev[deriv_func] * cs_d_fwd[deriv_func]

                    # fwd_sum = np.sum(func_seeds_fwd[deriv_func])

                    # # print(geom_seeds_rev)
                    tol = 1e-13
                    # print(f"{deriv_func} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(fwd_sum) < tol or np.abs(rev_sum) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            atol=1e-14,
                            err_msg=f"deriv_func {deriv_func} w.r.t. {surf_key}:{geom_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            rtol=1e-12,
                            err_msg=f"deriv_func {deriv_func} w.r.t. {surf_key}:{geom_key}",
                        )

    def test_fwd_gamma_d(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        num_consurf = self.ovl_solver.get_num_control_surfs()
        gamma_d_seeds = np.random.rand(num_consurf, num_gamma)

        cs_d = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds)[2]
        cs_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds, mode="FD", step=1e-7)[2]

        for deriv_func in cs_d:
            sens_label = f"{deriv_func} wrt gamma_d"
            # print(sens_label, cs_d[deriv_func], cs_d_fd[deriv_func])
            np.testing.assert_allclose(
                cs_d[deriv_func],
                cs_d_fd[deriv_func],
                rtol=1e-3,
                err_msg=sens_label,
            )

    def test_rev_gamma_d(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        num_consurf = self.ovl_solver.get_num_control_surfs()
        gamma_d_seeds_fwd = np.random.rand(num_consurf, num_gamma)

        cs_d_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds_fwd)[2]
        self.ovl_solver.clear_ad_seeds_fast()

        for deriv_func in cs_d_fwd:
                cs_d_rev = {deriv_func: 1.0}

                gamma_d_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_d_rev)[3]

                rev_sum = np.sum(gamma_d_seeds_rev * gamma_d_seeds_fwd)

                fwd_sum = np.sum(cs_d_fwd[deriv_func])

                # print("fwd_sum", fwd_sum, "rev_sum", rev_sum)
                np.testing.assert_allclose(
                    fwd_sum,
                    rev_sum,
                    atol=1e-14,
                    err_msg=f"deriv_func {deriv_func} w.r.t. gamma",
                )


if __name__ == "__main__":
    unittest.main()