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


class TestStabDerivs(unittest.TestCase):
    # TODO: beta derivatives likely wrong

    def setUp(self):
        # self.ovl_solver = OVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.ovl_solver = OVLSolver(geo_file="aircraft_L1.avl", debug=False)
        # self.ovl_solver = OVLSolver(geo_file="rect.avl")
        
        # HACK: we have to use alpha/beta=0 here so the stability and body axis are the same
        self.ovl_solver.set_constraint("alpha", 0.0)
        self.ovl_solver.set_constraint("beta", 0.0)
        
        self.ovl_solver.execute_run()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id()} Memory usage: {mb_memory:.2f} MB")

    def finite_dif(self, con_list, geom_seeds, step=1e-7):
        con_seeds = {}

        for con in con_list:
            con_seeds[con] = 1.0

        self.ovl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
        self.ovl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)

        self.ovl_solver.avl.update_surfaces()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.exec_rhs()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.velsum()
        self.ovl_solver.avl.aero()
        # self.ovl_solver.execute_run()
        coef_data_peturb = self.ovl_solver.get_total_forces()
        consurf_derivs_peturb = self.ovl_solver.get_control_body_axis_derivs()

        self.ovl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-step)
        self.ovl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-step)

        self.ovl_solver.execute_run()

        coef_data = self.ovl_solver.get_total_forces()
        consurf_derivs = self.ovl_solver.get_control_body_axis_derivs()

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

        return func_seeds, consurf_derivs_seeds

    def test_deriv_values(self):
        # compare the analytical gradients with finite difference for each constraint and function
        import pprint
        base_data = self.ovl_solver.get_total_forces()
        # pprint.pprint(base_data)
        
        # for key in base_data:
        #     print(f"{key:5} {base_data[key]}")
            
        
        
        body_axis_derivs = self.ovl_solver.get_body_axis_derivs()
        # pprint.pprint(body_axis_derivs)

        con_keys =  ["roll rate", "pitch rate", "yaw rate"]
        func_keys = ["CX body axis", "CY body axis", "CZ body axis", "Cl", "Cm", "Cn"]
        # con_keys =  ["roll rate"]
        # func_keys = ["CZ body axis"]
    
        func_to_key = {
            "Cl" : "Cl",
            "Cm" : "Cm",
            "Cn" : "Cn",
            "CX body axis" : "CX",
            "CY body axis" : "CY",
            "CZ body axis" : "CZ",
        }
        
        con_to_var = {
            "roll rate" : "p" ,
            "pitch rate" : "q" ,
            "yaw rate" : "r" ,
        }
        
        for con_key in con_keys:
            h = 1e-10
            val = self.ovl_solver.get_constraint(con_key)
            self.ovl_solver.set_constraint(con_key, val + h)
            self.ovl_solver.execute_run()
            perb_data = self.ovl_solver.get_total_forces()
            self.ovl_solver.set_constraint(con_key, val)
            
            for func_key in func_keys:
                key = self.ovl_solver._get_deriv_key(con_to_var[con_key], func_to_key[func_key])
                ad_dot = body_axis_derivs[key] 
                
                fd_dot = (perb_data[func_key] - base_data[func_key]) / h 
                
                if con_key in ['alpha', 'beta']:
                   fd_dot *= 180/np.pi 
                   # convert to radians from degrees!
                
                rel_err = np.abs((ad_dot - fd_dot) / (fd_dot + 1e-20))
                print(f"{key:5}  | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-8
                if np.abs(ad_dot) < tol or np.abs(fd_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        atol=tol,
                        err_msg=f"func_key {key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        rtol=5e-5,
                        err_msg=f"func_key {key}",
                    )   
            
            


if __name__ == "__main__":
    unittest.main()