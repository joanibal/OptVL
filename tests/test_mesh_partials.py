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
import pickle


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, '..', 'geom_files')
rect_file = os.path.join(geom_dir, 'rect_with_body.pkl')
with open(rect_file, 'rb') as f:
    input_dict = pickle.load(f)


class TestFunctionPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=rect_file)
        self.ovl_solver = OVLSolver(input_dict=input_dict)
        self.ovl_solver.set_variable("alpha", 25.0)
        self.ovl_solver.set_variable("beta", 5.0)
        self.ovl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id():80} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_mesh(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                func_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}
                )[0]

                func_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}, mode="FD", step=1e-7
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
                            err_msg=f"func_key {func_key} w.r.t. {mesh_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            func_seeds[func_key],
                            func_seeds_FD[func_key],
                            rtol=5e-4,
                            err_msg=f"func_key {func_key} w.r.t. {mesh_key}",
                        )

    def test_rev_mesh(self):
        np.random.seed(111)

        sens_dict_rev = {}
        for func_key in self.ovl_solver.case_var_to_fort_var:
            sens_dict_rev[func_key] = self.ovl_solver._execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})[2]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                func_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}
                )[0]

                for func_key in func_seeds_fwd:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(sens_dict_rev[func_key][surf_key][mesh_key] * mesh_seeds)
                    fwd_sum = np.sum(func_seeds_fwd[func_key])

                    # # print(mesh_seeds_rev)
                    tol = 1e-13
                    # print(f"{func_key} wrt {surf_key}:{mesh_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(fwd_sum) < tol or np.abs(rev_sum) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            atol=1e-14,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{mesh_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            rtol=1e-12,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{mesh_key}",
                        )

class TestResidualPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file)
        self.ovl_solver = OVLSolver(input_dict=input_dict)
        self.ovl_solver.set_variable("alpha", 25.0)
        self.ovl_solver.set_variable("beta", 5.0)
        self.ovl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id():80} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_mesh(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}
                )[1]

                res_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}, mode="FD", step=1e-8
                )[1]

                abs_error = np.abs(res_seeds - res_seeds_FD)
                rel_error = (res_seeds - res_seeds_FD) / (res_seeds + 1e-15)
                idx_max_rel_error = np.argmax(np.abs(rel_error))
                idx_max_abs_error = np.argmax(np.abs(abs_error))

                # print(
                #     f"{surf_key:10} {mesh_key:10} AD:{np.linalg.norm(res_seeds): .5e} FD:{np.linalg.norm(res_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                # )
                np.testing.assert_allclose(
                    res_seeds,
                    res_seeds_FD,
                    atol=3e-5,
                    err_msg=f"func_key res w.r.t. {surf_key}:{mesh_key}",
                )

    def test_rev_mesh(self):
        num_res = self.ovl_solver.get_mesh_size()
        res_seeds_rev = np.random.seed(111)
        res_seeds_rev = np.random.rand(num_res)
        mesh_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)[2]

        self.ovl_solver.clear_ad_seeds_fast()
        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                res_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}
                )[1]

                # do dot product
                res_sum = np.sum(res_seeds_rev * res_seeds)
                geom_sum = np.sum(mesh_seeds_rev[surf_key][mesh_key] * mesh_seeds)

                # print(f"res wrt {surf_key}:{mesh_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    geom_sum,
                    atol=1e-14,
                    err_msg=f"func_key res w.r.t. {surf_key}:{mesh_key}",
                )

class TestResidualUPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file)
        self.ovl_solver = OVLSolver(input_dict=input_dict)
        self.ovl_solver.set_variable("alpha", 25.0)
        self.ovl_solver.set_variable("beta", 5.0)
        self.ovl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_mesh(self):
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                res_u_seeds = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}
                )[6]

                res_u_seeds_FD = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}, mode="FD", step=1e-8
                )[6]

                abs_error = np.abs(res_u_seeds.flatten() - res_u_seeds_FD.flatten())
                rel_error = np.abs((res_u_seeds.flatten() - res_u_seeds_FD.flatten()) / (res_u_seeds.flatten() + 1e-15))

                idx_max_rel_error = np.argmax(rel_error)
                idx_max_abs_error = np.argmax(abs_error)
                # print(
                #     f"{surf_key:10} {mesh_key:10} AD:{np.linalg.norm(res_u_seeds): .5e} FD:{np.linalg.norm(res_u_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                # )
                np.testing.assert_allclose(
                    res_u_seeds,
                    res_u_seeds_FD,
                    atol=1e-4,
                    err_msg=f" res_u w.r.t. {surf_key}:{mesh_key}",
                )

    def test_rev_mesh(self):
        np.random.seed(111)
        num_gamma = self.ovl_solver.get_mesh_size()
        res_u_seeds_rev = np.random.rand(self.ovl_solver.NUMAX, num_gamma)

        self.ovl_solver.clear_ad_seeds_fast()

        mesh_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(res_u_seeds=res_u_seeds_rev)[2]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                res_u_seeds_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds}}
                )[6]

                # do dot product
                res_sum = np.sum(res_u_seeds_rev * res_u_seeds_fwd)
                mesh_sum = np.sum(mesh_seeds_rev[surf_key][mesh_key] * mesh_seeds)

                # print(f"res wrt {surf_key}:{mesh_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    mesh_sum,
                    atol=1e-14,
                    err_msg=f"res_u w.r.t. {surf_key}:{mesh_key}",
                )
        self.ovl_solver.clear_ad_seeds_fast()

class TestStabDerivDerivsPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        # self.ovl_solver = OVLSolver(geo_file=geom_file)
        self.ovl_solver = OVLSolver(input_dict=input_dict)
        # self.ovl_solver = OVLSolver(geo_file="geom_files/rect.avl")
        self.ovl_solver.set_variable("alpha", 45.0)
        self.ovl_solver.set_variable("beta", 45.0)
        self.ovl_solver.execute_run()
        self.ovl_solver.clear_ad_seeds_fast()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_mesh(self):
        # this one is broken start here
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds = np.random.rand(*arr.shape)

                sd_d = self.ovl_solver._execute_jac_vec_prod_fwd(mesh_seeds={surf_key: {mesh_key: mesh_seeds}})[3]

                sd_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    mesh_seeds={surf_key: {mesh_key: mesh_seeds}}, mode="FD", step=5e-8
                )[3]

                for deriv_func in sd_d:
                    sens_label = f"{deriv_func} wrt {surf_key}:{mesh_key:5}"

                    # print(f"{sens_label} AD:{sd_d[deriv_func]} FD:{sd_d_fd[deriv_func]}")
                    # quit()
                    tol = 1e-10
                    # print(f"{deriv_func} wrt {surf_key}:{mesh_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(sd_d[deriv_func]) < tol or np.abs(sd_d_fd[deriv_func]) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        # this is basiccally saying if one is less than 1e-10 the other must be less than 5e-7
                        np.testing.assert_allclose(
                            sd_d[deriv_func],
                            sd_d_fd[deriv_func],
                            atol=5e-7,
                            err_msg=sens_label,
                        )
                    else:
                        np.testing.assert_allclose(
                            sd_d[deriv_func],
                            sd_d_fd[deriv_func],
                            rtol=5e-3,
                            err_msg=sens_label,
                        )

    def test_rev_mesh(self):
        np.random.seed(111)
        sd_d_rev = {}
        for deriv_func in self.ovl_solver.case_stab_derivs_to_fort_var:
            sd_d_rev[deriv_func] = np.random.rand(1)[0]

        mesh_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(stab_derivs_seeds=sd_d_rev)[2]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_mesh_to_fort_var:
            for mesh_key in self.ovl_solver.surf_mesh_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                mesh_seeds_fwd = np.random.rand(*arr.shape)

                sd_d_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, mesh_seeds={surf_key: {mesh_key: mesh_seeds_fwd}}
                )[3]

                for deriv_func in self.ovl_solver.case_stab_derivs_to_fort_var:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(mesh_seeds_rev[surf_key][mesh_key] * mesh_seeds_fwd)

                    fwd_sum = 0.0
                    for deriv_func in sd_d_fwd:
                        fwd_sum += sd_d_rev[deriv_func] * sd_d_fwd[deriv_func]

                    # # print(mesh_seeds_rev)
                    tol = 1e-13
                    # print(f"{deriv_func} wrt {surf_key}:{mesh_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(fwd_sum) < tol or np.abs(rev_sum) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            atol=1e-14,
                            err_msg=f"deriv_func {deriv_func} w.r.t. {surf_key}:{mesh_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            rtol=1e-12,
                            err_msg=f"deriv_func {deriv_func} w.r.t. {surf_key}:{mesh_key}",
                        )


if __name__ == "__main__":
    unittest.main()
