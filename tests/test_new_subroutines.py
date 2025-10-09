# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver
import copy

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


class TestNewSubroutines(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl", debug=False)
        self.ovl_solver.set_constraint("alpha", 25.0)

    def test_residual(self):
        self.ovl_solver.avl.get_res()
        res = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RES)
        rhs = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RHS)
        np.testing.assert_allclose(
            res,
            -1 * rhs,
            atol=1e-15,
        )
        res_u = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RES_U)
        rhs_u = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RHS_U)

        np.testing.assert_allclose(
            res_u,
            -1 * rhs_u,
            atol=1e-15,
        )

        self.ovl_solver.avl.exec_rhs()
        self.ovl_solver.avl.get_res()
        res = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RES)
        res_d = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RES_D)
        res_u = copy.deepcopy(self.ovl_solver.avl.VRTX_R.RES_U)

        np.testing.assert_allclose(
            res,
            np.zeros_like(res),
            atol=5e-14,
        )

        np.testing.assert_allclose(
            res_d,
            np.zeros_like(res_d),
            atol=5e-14,
        )

        np.testing.assert_allclose(
            res_u,
            np.zeros_like(res_u),
            atol=5e-14,
        )

    def test_new_solve(self):
        self.ovl_solver.set_constraint("Elevator", 10.00)
        self.ovl_solver.set_constraint("alpha", 10.00)
        self.ovl_solver.set_constraint("beta", 10.00)

        self.ovl_solver.avl.exec_rhs()

        self.ovl_solver.avl.velsum()
        self.ovl_solver.avl.aero()
        wv_new = copy.deepcopy(self.ovl_solver.avl.SOLV_R.WV)
        gam_new = copy.deepcopy(self.ovl_solver.avl.VRTX_R.GAM)
        gam_u_new = copy.deepcopy(self.ovl_solver.avl.VRTX_R.GAM_U)
        coef_data_new = self.ovl_solver.get_total_forces()

        coef_derivs_new = self.ovl_solver.get_control_stab_derivs()

        self.ovl_solver.execute_run()
        gam = copy.deepcopy(self.ovl_solver.avl.VRTX_R.GAM)
        wv = copy.deepcopy(self.ovl_solver.avl.SOLV_R.WV)
        gam_u = copy.deepcopy(self.ovl_solver.avl.VRTX_R.GAM_U)
        coef_data = self.ovl_solver.get_total_forces()
        coef_derivs = self.ovl_solver.get_control_stab_derivs()

        np.testing.assert_allclose(
            wv,
            wv_new,
            atol=1e-15,
        )
        np.testing.assert_allclose(
            gam,
            gam_new,
            atol=1e-15,
        )
        np.testing.assert_allclose(
            gam_u,
            gam_u_new,
            atol=1e-15,
        )
        for func_key in coef_data:
            np.testing.assert_allclose(
                coef_data[func_key],
                coef_data_new[func_key],
                atol=1e-14,
                err_msg=f"func_key {func_key}",
            )

        for func_key in coef_derivs:
            # for consurf_key in coef_derivs[func_key]:
            np.testing.assert_allclose(
                coef_derivs[func_key],
                coef_derivs_new[func_key],
                err_msg=f"deriv of func_key {func_key}",
                atol=1e-14,
            )


if __name__ == "__main__":
    unittest.main()
