# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver, OVLGroup

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np
import openmdao.api as om
import warnings

# Set DeprecationWarning to be treated as an error
warnings.simplefilter("error", DeprecationWarning)

# Optionally, you can also set NumPy options to raise errors for specific conditions
np.seterr(all="raise")

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")


class TestOMWrapper(unittest.TestCase):
    def setUp(self):
        self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)

        model = om.Group()
        model.add_subsystem(
            "ovlsolver",
            OVLGroup(
                geom_file=geom_file,
                mass_file=mass_file,
                output_stability_derivs=True,
                input_param_vals=True,
                input_ref_vals=True,
                input_airfoil_geom=True,
            ),
        )

        self.prob = om.Problem(model)

    def test_aero_coef(self):
        self.ovl_solver.execute_run()
        run_data = self.ovl_solver.get_total_forces()

        prob = self.prob
        prob.setup(mode="rev")
        prob.run_model()

        for func in run_data:
            om_val = prob.get_val(f"ovlsolver.{func}")
            assert om_val == run_data[func]

    def test_surface_param_setting(self):
        prob = self.prob
        prob.setup(mode="rev")

        np.random.seed(111)

        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            for geom_key in self.ovl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.ovl_solver.get_surface_param(surf_key, geom_key)
                arr += np.random.rand(*arr.shape) * 0.1
                # print(f'setting {surf_key}:{geom_key} to {arr}')
                # #set surface data
                self.ovl_solver.set_surface_param(surf_key, geom_key, arr)
                self.ovl_solver.avl.update_surfaces()
                self.ovl_solver.execute_run()

                # set om surface data
                prob.set_val(f"ovlsolver.{surf_key}:{geom_key}", arr)
                prob.run_model()

                run_data = self.ovl_solver.get_total_forces()
                for func in run_data:
                    om_val = prob.get_val(f"ovlsolver.{func}")
                    assert om_val == run_data[func]

                stab_derivs = self.ovl_solver.get_stab_derivs()
                for func in stab_derivs:
                    om_val = prob.get_val(f"ovlsolver.{func}")
                    assert om_val == stab_derivs[func]


    def test_param_setting(self):
        prob = self.prob
        prob.setup(mode="rev")

        np.random.seed(111)

        for param in ["CD0", "Mach", "X cg", "Y cg", "Z cg"]:
            param_val = self.ovl_solver.get_parameter(param)
            param_val += np.random.rand() * 0.1
            # #set surface data
            self.ovl_solver.set_parameter(param, param_val)
            self.ovl_solver.execute_run()

            # set om surface data
            prob.set_val(f"ovlsolver.{param}", param_val)
            prob.run_model()

            run_data = self.ovl_solver.get_total_forces()
            for func in run_data:
                om_val = prob.get_val(f"ovlsolver.{func}")
                assert om_val == run_data[func]

            stab_derivs = self.ovl_solver.get_stab_derivs()
            for func in stab_derivs:
                om_val = prob.get_val(f"ovlsolver.{func}")
                print(func, om_val, stab_derivs[func])
                assert om_val == stab_derivs[func]


    def test_CL_solve(self):
        prob = self.prob
        cl_star = 1.5
        prob.model.add_design_var("ovlsolver.alpha", lower=-10, upper=10)
        prob.model.add_constraint("ovlsolver.CL", equals=cl_star)
        prob.model.add_objective("ovlsolver.CD", scaler=1e3)
        prob.setup(mode="rev")
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["debug_print"] = ["desvars", "ln_cons", "nl_cons", "objs"]
        prob.driver.options["tol"] = 1e-6
        prob.driver.options["disp"] = True

        prob.setup(mode="rev")
        prob.run_driver()
        om.n2(prob, show_browser=False, outfile="vlm_opt.html")

        om_val = prob.get_val("ovlsolver.alpha")

        self.ovl_solver.set_trim_condition("CL", cl_star)
        self.ovl_solver.execute_run()
        alpha = self.ovl_solver.get_parameter("alpha")

        np.testing.assert_allclose(
            om_val,
            alpha,
            rtol=1e-5,
            err_msg="solved alpha",
        )

    def test_CM_solve(self):
        prob = self.prob
        prob.model.add_design_var("ovlsolver.alpha", lower=-10, upper=10)
        prob.model.add_constraint("ovlsolver.CM", equals=0.0, scaler=1e3)
        prob.model.add_objective("ovlsolver.CD", scaler=1e3)
        prob.setup(mode="rev")
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options["optimizer"] = "SLSQP"
        prob.driver.options["debug_print"] = ["desvars", "ln_cons", "nl_cons", "objs"]
        prob.driver.options["tol"] = 1e-6
        prob.driver.options["disp"] = True

        prob.setup(mode="rev")
        prob.run_driver()
        om.n2(prob, show_browser=False, outfile="vlm_opt.html")

        om_val = prob.get_val("ovlsolver.alpha")

        self.ovl_solver.set_constraint("alpha", 0.00, con_var="Cm pitch moment")
        self.ovl_solver.execute_run()
        alpha = self.ovl_solver.get_parameter("alpha")

        np.testing.assert_allclose(
            om_val,
            alpha,
            rtol=1e-5,
            err_msg="solved alpha",
        )

    def test_OM_total_derivs(self):
        prob = self.prob
        cl_star = 1.5
        dcl_dalpha_star = -0.1
        prob.model.add_design_var("ovlsolver.Wing:xles")
        prob.model.add_design_var("ovlsolver.Wing:yles")
        prob.model.add_design_var("ovlsolver.Wing:zles")
        prob.model.add_design_var("ovlsolver.Wing:chords")
        prob.model.add_design_var("ovlsolver.Wing:aincs")
        prob.model.add_design_var("ovlsolver.Elevator", lower=-10, upper=10)
        prob.model.add_design_var("ovlsolver.alpha", lower=-10, upper=10)
        prob.model.add_design_var("ovlsolver.Sref")
        prob.model.add_design_var("ovlsolver.Mach")
        prob.model.add_design_var("ovlsolver.X cg")
        prob.model.add_constraint("ovlsolver.CL", equals=cl_star)
        prob.model.add_constraint("ovlsolver.dCL/dalpha", equals=-dcl_dalpha_star)
        prob.model.add_objective("ovlsolver.CD", scaler=1e3)
        prob.model.add_objective("ovlsolver.CM", scaler=1e3)
        prob.setup(mode="rev")
        prob.run_model()
        om.n2(prob, show_browser=False, outfile="vlm_opt.html")

        deriv_err = prob.check_totals()
        rtol = 5e-4
        for key, data in deriv_err.items():
            np.testing.assert_allclose(
                data["J_fd"],
                data["J_rev"],
                rtol=rtol,
                err_msg=f"deriv of {key[0]} wrt {key[1]} does not agree with FD to rtol={rtol}",
            )


if __name__ == "__main__":
    unittest.main()
