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


class TestBodyAxisDerivDerivsPartials(unittest.TestCase):
    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl")
        # self.ovl_solver = OVLSolver(geo_file="rect.avl")
        self.ovl_solver.set_variable("alpha", 45.0)
        self.ovl_solver.set_variable("beta", 45.0)
        self.ovl_solver.execute_run()
        self.ovl_solver.clear_ad_seeds_fast()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def test_fwd_aero_constraint(self):
        for con_key in self.ovl_solver.con_var_to_fort_var:
            bd_d = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[4]

            bd_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, mode="FD", step=1e-6)[4]

            for deriv_func in bd_d:
                sens_label = f"{deriv_func} wrt {con_key}"
                # print(sens_label, bd_d[deriv_func], bd_d_fd[deriv_func])
                np.testing.assert_allclose(
                    bd_d[deriv_func],
                    bd_d_fd[deriv_func],
                    rtol=1e-4,
                    err_msg=sens_label,
                )

    def test_rev_aero_constraint(self):
        
        body_axis_deriv_seeds_rev = {}
        for deriv_func, var_dict in self.ovl_solver.case_body_derivs_to_fort_var.items():
            body_axis_deriv_seeds_rev[deriv_func] = np.random.rand(1)[0]

        con_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(body_axis_derivs_seeds=body_axis_deriv_seeds_rev)[0]

        self.ovl_solver.clear_ad_seeds_fast()

        for con_key in self.ovl_solver.con_var_to_fort_var:
            body_axis_deriv_seeds_fwd= self.ovl_solver._execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})[4]

            body_axis_deriv_sum = 0.0
            for deriv_func in body_axis_deriv_seeds_rev:
                    body_axis_deriv_sum += body_axis_deriv_seeds_rev[deriv_func] * body_axis_deriv_seeds_fwd[deriv_func]

            # do dot product
            con_sum = np.sum(con_seeds_rev[con_key])

            # print(f"cs_dervs wrt {con_key}", "rev", con_sum, "fwd", body_axis_deriv_sum)

            np.testing.assert_allclose(
                con_sum,
                body_axis_deriv_sum,
                atol=1e-14,
                err_msg=f"cs_dervs wrt {con_key}",
            )

    def test_fwd_geom(self):
        # this one is broken start here
        np.random.seed(111)
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                bd_d = self.ovl_solver._execute_jac_vec_prod_fwd(geom_seeds={surf_key: {geom_key: geom_seeds}})[4]

                bd_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-9
                )[4]

                for deriv_func in bd_d:
                    sens_label = f"{deriv_func} wrt {surf_key}:{geom_key:5}"

                    # print(f"{sens_label} AD:{bd_d[deriv_func]} FD:{bd_d_fd[deriv_func]}")
                    # quit()
                    tol = 5e-7
                    # print(f"{deriv_func} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(bd_d[deriv_func]) < tol or np.abs(bd_d_fd[deriv_func]) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            bd_d[deriv_func],
                            bd_d_fd[deriv_func],
                            atol=tol,
                            err_msg=sens_label,
                        )
                    else:
                        np.testing.assert_allclose(
                            bd_d[deriv_func],
                            bd_d_fd[deriv_func],
                            rtol=5e-3,
                            err_msg=sens_label,
                        )

    def test_rev_geom(self):
        np.random.seed(111)
        bd_d_rev = {}
        for deriv_func in self.ovl_solver.case_body_derivs_to_fort_var:
            bd_d_rev[deriv_func] = np.random.rand(1)[0]


        geom_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(body_axis_derivs_seeds=bd_d_rev)[1]
        self.ovl_solver.clear_ad_seeds_fast()

        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds_fwd = np.random.rand(*arr.shape)

                bd_d_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds_fwd}}
                )[4]

                for deriv_func in self.ovl_solver.case_body_derivs_to_fort_var:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds_fwd)

                    fwd_sum = 0.0
                    for deriv_func in bd_d_fwd:
                        fwd_sum += bd_d_rev[deriv_func] * bd_d_fwd[deriv_func]


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
        num_gamma = self.ovl_solver.get_mesh_size()
        gamma_u_seeds = np.random.rand(self.ovl_solver.NUMAX, num_gamma)

        bd_d = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds)[4]
        bd_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds, mode="FD", step=1e-7)[4]

        for deriv_func in bd_d:
            sens_label = f"{deriv_func} wrt gamma_u"
            np.testing.assert_allclose(
                bd_d[deriv_func],
                bd_d_fd[deriv_func],
                rtol=1e-4,
                err_msg=sens_label,
            )

    def test_rev_gamma_u(self):
        num_gamma = self.ovl_solver.get_mesh_size()
        num_consurf = self.ovl_solver.get_num_control_surfs()
        gamma_u_seeds_fwd = np.random.rand(num_consurf, num_gamma)

        bd_d_fwd = self.ovl_solver._execute_jac_vec_prod_fwd(gamma_u_seeds=gamma_u_seeds_fwd)[4]
        self.ovl_solver.clear_ad_seeds_fast()


        for deriv_func in bd_d_fwd:
            # for var_key in bd_d_fwd[deriv_func]:
            bd_d_rev = {deriv_func: 1.0}

            gamma_u_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(body_axis_derivs_seeds=bd_d_rev)[4]

            rev_sum = np.sum(gamma_u_seeds_rev * gamma_u_seeds_fwd)

            fwd_sum = np.sum(bd_d_fwd[deriv_func])

            # print("fwd_sum", fwd_sum, "rev_sum", rev_sum)
            np.testing.assert_allclose(
                fwd_sum,
                rev_sum,
                atol=1e-14,
                err_msg=f"deriv_func {deriv_func} w.r.t. gamma",
            )


    def test_fwd_ref(self):
        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            bd_d = self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[4]

            bd_d_fd = self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0}, mode="FD", step=1e-6)[4]

            for deriv_func in bd_d:
                sens_label = f"{deriv_func} wrt {ref_key}"
                # print(sens_label, bd_d[deriv_func], bd_d_fd[deriv_func])
                
                tol = 1e-8
                if np.abs( bd_d[deriv_func]) < tol or np.abs(bd_d_fd[deriv_func]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        bd_d[deriv_func],
                        bd_d_fd[deriv_func],
                        atol=1e-5,
                        err_msg=sens_label,
                    )
                else:
                    np.testing.assert_allclose(
                        bd_d[deriv_func],
                        bd_d_fd[deriv_func],
                        rtol=1e-5,
                        err_msg=sens_label,
                    )
                        


    def test_rev_ref(self):
        
        body_axis_deriv_seeds_rev = {}
        for deriv_func, var_dict in self.ovl_solver.case_body_derivs_to_fort_var.items():
            body_axis_deriv_seeds_rev[deriv_func] = np.random.rand(1)[0]

        ref_seeds_rev = self.ovl_solver._execute_jac_vec_prod_rev(body_axis_derivs_seeds=body_axis_deriv_seeds_rev)[6]

        self.ovl_solver.clear_ad_seeds_fast()

        for ref_key in self.ovl_solver.ref_var_to_fort_var:
            body_axis_deriv_seeds_fwd= self.ovl_solver._execute_jac_vec_prod_fwd(ref_seeds={ref_key: 1.0})[4]

            body_axis_deriv_sum = 0.0
            for deriv_func in body_axis_deriv_seeds_fwd:
                    body_axis_deriv_sum += body_axis_deriv_seeds_rev[deriv_func] * body_axis_deriv_seeds_fwd[deriv_func]

            # do dot product
            ref_sum = np.sum(ref_seeds_rev[ref_key])

            # print(f"cs_dervs wrt {ref_key}", "rev", ref_sum, "fwd", body_axis_deriv_sum)

            np.testing.assert_allclose(
                ref_sum,
                body_axis_deriv_sum,
                atol=1e-14,
                err_msg=f"cs_dervs wrt {ref_key}",
            )

                
if __name__ == "__main__":
    unittest.main()
