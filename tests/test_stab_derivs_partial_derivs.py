# =============================================================================
# Extension modules
# =============================================================================
from optvl import AVLSolver
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

class TestResidualUPartials(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.set_constraint("alpha", 25.0)
        self.avl_solver.set_constraint("beta", 5.0)
        self.avl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")


    def test_fwd_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            res_u_seeds = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[5]

            res_u_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, mode="FD", step=1e-5
            )[5]

            np.testing.assert_allclose(
                res_u_seeds,
                res_u_seeds_FD,
                rtol=1e-5,
            )

    def test_rev_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
        # for con_key in ["beta", "beta"]:
            res_u_seeds = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[5]

            num_gamma = self.avl_solver.get_mesh_size()
            np.random.seed(111)
            res_u_seeds_rev = np.random.rand(self.avl_solver.NUMAX, num_gamma)

            self.avl_solver.clear_ad_seeds_fast()

            con_seeds = self.avl_solver.execute_jac_vec_prod_rev(res_u_seeds=res_u_seeds_rev)[0]
            self.avl_solver.clear_ad_seeds_fast()

            # do dot product
            res_sum = np.sum(res_u_seeds_rev * res_u_seeds)

            # print(f"res wrt {con_key}", "fwd", res_sum, "rev", con_seeds[con_key])

            np.testing.assert_allclose(
                res_sum,
                con_seeds[con_key],
                atol=1e-14,
                err_msg=f"deriv_func res w.r.t. {con_key}",
            )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                res_u_seeds = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[5]

                res_u_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )[5]

                abs_error = np.abs(res_u_seeds.flatten() - res_u_seeds_FD.flatten())
                rel_error = np.abs((res_u_seeds.flatten() - res_u_seeds_FD.flatten()) / (res_u_seeds.flatten() + 1e-15))

                idx_max_rel_error = np.argmax(rel_error)
                idx_max_abs_error = np.argmax(abs_error)
                # print(
                #     f"{surf_key:10} {geom_key:10} AD:{np.linalg.norm(res_u_seeds): .5e} FD:{np.linalg.norm(res_u_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                # )
                np.testing.assert_allclose(
                    res_u_seeds,
                    res_u_seeds_FD,
                    atol=1e-4,
                    err_msg=f" res_u w.r.t. {surf_key}:{geom_key}",
                )

    def test_rev_geom(self):
        np.random.seed(111)
        num_gamma = self.avl_solver.get_mesh_size()
        res_u_seeds_rev = np.random.rand(self.avl_solver.NUMAX, num_gamma)

        self.avl_solver.clear_ad_seeds_fast()

        geom_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(res_u_seeds=res_u_seeds_rev)[1]
        self.avl_solver.clear_ad_seeds_fast()

        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                res_u_seeds_fwd = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )[5]

                # do dot product
                res_sum = np.sum(res_u_seeds_rev * res_u_seeds_fwd)
                geom_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds)

                # print(f"res wrt {surf_key}:{geom_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    geom_sum,
                    atol=1e-14,
                    err_msg=f"res_u w.r.t. {surf_key}:{geom_key}",
                )
        self.avl_solver.clear_ad_seeds_fast()

    def test_fwd_gamma_u(self):
        # stoped here for the night
        num_gamma = self.avl_solver.get_mesh_size()
        gamma_u_seeds = np.random.rand(self.avl_solver.NUMAX, num_gamma)
        # gamma_u_seeds = np.array([[1],[0]])

        res_u_seeds = self.avl_solver.execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds)[5]

        res_u_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
            gamma_u_seeds=gamma_u_seeds, mode="FD", step=1e-0
        )[5]

        np.testing.assert_allclose(
            res_u_seeds,
            res_u_seeds_FD,
            rtol=1e-5,
        )

    def test_rev_gamma_u(self):
        num_gamma = self.avl_solver.get_mesh_size()
        gamma_u_seeds_fwd = np.random.rand(self.avl_solver.NUMAX, num_gamma)

        res_u_seeds_rev = np.random.rand(self.avl_solver.NUMAX, num_gamma)

        res_u_seeds_fwd = self.avl_solver.execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds_fwd)[5]
        self.avl_solver.clear_ad_seeds_fast()

        gamma_u_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(res_u_seeds=res_u_seeds_rev)[4]

        gamma_sum = np.sum(gamma_u_seeds_rev * gamma_u_seeds_fwd)
        res_sum = np.sum(res_u_seeds_rev * res_u_seeds_fwd)

        # print("fwd_sum", gamma_sum, "rev_sum", res_sum)
        np.testing.assert_allclose(
            gamma_sum,
            res_sum,
            atol=1e-14,
            err_msg=f"res w.r.t. gamma",
        )
    
    def test_fwd_ref(self):
        for ref_key in self.avl_solver.ref_var_to_fort_var:
            res_u_seeds = self.avl_solver.execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[5]

            res_u_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
                ref_seeds={ref_key: 1.0}, mode="FD", step=1e-5
            )[5]
            
            print(res_u_seeds, res_u_seeds_FD)

            np.testing.assert_allclose(
                res_u_seeds,
                res_u_seeds_FD,
                rtol=1e-5,
            )



class TestStabDerivDerivsPartials(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.set_constraint("alpha", 45.0)
        self.avl_solver.set_constraint("beta", 45.0)
        self.avl_solver.execute_run()
        self.avl_solver.clear_ad_seeds_fast()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            sd_d = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[3]

            sd_d_fd = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, mode="FD", step=1e-6)[3]

            for deriv_func in sd_d:
                sens_label = f"{deriv_func} wrt {con_key}"
                # print(sens_label, sd_d[deriv_func][cs_key], sd_d_fd[deriv_func][cs_key])
                np.testing.assert_allclose(
                    sd_d[deriv_func],
                    sd_d_fd[deriv_func],
                    rtol=1e-4,
                    err_msg=sens_label,
                )

    def test_rev_aero_constraint(self):
        
        stab_deriv_seeds_rev = {}
        for deriv_func, var_dict in self.avl_solver.case_stab_derivs_to_fort_var.items():
            stab_deriv_seeds_rev[deriv_func] = np.random.rand(1)[0]

        con_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(stab_derivs_seeds=stab_deriv_seeds_rev)[0]

        self.avl_solver.clear_ad_seeds_fast()

        for con_key in self.avl_solver.con_var_to_fort_var:
            stab_deriv_seeds_fwd= self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[3]

            stab_deriv_sum = 0.0
            for deriv_func in stab_deriv_seeds_fwd:
                    stab_deriv_sum += stab_deriv_seeds_rev[deriv_func] * stab_deriv_seeds_fwd[deriv_func]

            # do dot product
            con_sum = np.sum(con_seeds_rev[con_key])

            # print(f"cs_dervs wrt {con_key}", "rev", con_sum, "fwd", stab_deriv_sum)

            np.testing.assert_allclose(
                con_sum,
                stab_deriv_sum,
                atol=1e-14,
                err_msg=f"cs_dervs wrt {con_key}",
            )

    def test_fwd_geom(self):
        # this one is broken start here
        np.random.seed(111)
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                sd_d = self.avl_solver.execute_jac_vec_prod_fwd(geom_seeds={surf_key: {geom_key: geom_seeds}})[3]

                sd_d_fd = self.avl_solver.execute_jac_vec_prod_fwd(
                    geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=5e-8
                )[3]

                for deriv_func in sd_d:
                    sens_label = f"{deriv_func} wrt {surf_key}:{geom_key:5}"

                    # print(f"{sens_label} AD:{sd_d[deriv_func]} FD:{sd_d_fd[deriv_func]}")
                    # quit()
                    tol = 1e-10
                    # print(f"{deriv_func} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(sd_d[deriv_func]) < tol or np.abs(sd_d_fd[deriv_func]) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            sd_d[deriv_func],
                            sd_d_fd[deriv_func],
                            atol=1e-8,
                            err_msg=sens_label,
                        )
                    else:
                        np.testing.assert_allclose(
                            sd_d[deriv_func],
                            sd_d_fd[deriv_func],
                            rtol=1e-4,
                            err_msg=sens_label,
                        )

    def test_rev_geom(self):
        np.random.seed(111)
        sd_d_rev = {}
        for deriv_func in self.avl_solver.case_stab_derivs_to_fort_var:
            sd_d_rev[deriv_func] = np.random.rand(1)[0]


        geom_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(stab_derivs_seeds=sd_d_rev)[1]
        self.avl_solver.clear_ad_seeds_fast()

        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds_fwd = np.random.rand(*arr.shape)

                sd_d_fwd = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds_fwd}}
                )[3]

                for deriv_func in self.avl_solver.case_stab_derivs_to_fort_var:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds_fwd)

                    fwd_sum = 0.0
                    for deriv_func in sd_d_fwd:
                        fwd_sum += sd_d_rev[deriv_func] * sd_d_fwd[deriv_func]


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

    def test_fwd_gamma_u(self):
        num_gamma = self.avl_solver.get_mesh_size()
        gamma_u_seeds = np.random.rand(self.avl_solver.NUMAX, num_gamma)

        sd_d = self.avl_solver.execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds)[3]
        sd_d_fd = self.avl_solver.execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds, mode="FD", step=1e-7)[3]

        for deriv_func in sd_d:
            sens_label = f"{deriv_func} wrt gamma_u"
            np.testing.assert_allclose(
                sd_d[deriv_func],
                sd_d_fd[deriv_func],
                rtol=1e-6,
                err_msg=sens_label,
            )

    def test_rev_gamma_u(self):
        num_gamma = self.avl_solver.get_mesh_size()
        num_consurf = self.avl_solver.get_num_control_surfs()
        gamma_u_seeds_fwd = np.random.rand(num_consurf, num_gamma)

        sd_d_fwd = self.avl_solver.execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds_fwd)[3]
        self.avl_solver.clear_ad_seeds_fast()


        for deriv_func in sd_d_fwd:
            # for var_key in sd_d_fwd[deriv_func]:
            sd_d_rev = {deriv_func: 1.0}

            gamma_u_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(stab_derivs_seeds=sd_d_rev)[4]

            rev_sum = np.sum(gamma_u_seeds_rev * gamma_u_seeds_fwd)

            fwd_sum = np.sum(sd_d_fwd[deriv_func])

            # print("fwd_sum", fwd_sum, "rev_sum", rev_sum)
            np.testing.assert_allclose(
                fwd_sum,
                rev_sum,
                atol=1e-14,
                err_msg=f"deriv_func {deriv_func} w.r.t. gamma",
            )


    def test_fwd_ref(self):
        for ref_key in self.avl_solver.ref_var_to_fort_var:
            sd_d = self.avl_solver.execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[3]

            sd_d_fd = self.avl_solver.execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0}, mode="FD", step=1e-6)[3]

            for deriv_func in sd_d:
                sens_label = f"{deriv_func} wrt {ref_key}"
                print(sens_label, sd_d[deriv_func], sd_d_fd[deriv_func])
                np.testing.assert_allclose(
                    sd_d[deriv_func],
                    sd_d_fd[deriv_func],
                    rtol=1e-5,
                    err_msg=sens_label,
                )

    def test_rev_ref(self):
        
        stab_deriv_seeds_rev = {}
        for deriv_func, var_dict in self.avl_solver.case_stab_derivs_to_fort_var.items():
            stab_deriv_seeds_rev[deriv_func] = np.random.rand(1)[0]

        ref_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(stab_derivs_seeds=stab_deriv_seeds_rev)[6]

        self.avl_solver.clear_ad_seeds_fast()

        for ref_key in self.avl_solver.ref_var_to_fort_var:
            stab_deriv_seeds_fwd= self.avl_solver.execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[3]

            stab_deriv_sum = 0.0
            for deriv_func in stab_deriv_seeds_fwd:
                    stab_deriv_sum += stab_deriv_seeds_rev[deriv_func] * stab_deriv_seeds_fwd[deriv_func]

            # do dot product
            ref_sum = np.sum(ref_seeds_rev[ref_key])

            print(f"cs_dervs wrt {ref_key}", "rev", ref_sum, "fwd", stab_deriv_sum)

            np.testing.assert_allclose(
                ref_sum,
                stab_deriv_sum,
                atol=1e-14,
                err_msg=f"cs_dervs wrt {ref_key}",
            )

                
if __name__ == "__main__":
    unittest.main()
