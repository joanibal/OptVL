# =============================================================================
# Standard modules
# =============================================================================
import os
import copy
 
# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np
import sys
import psutil
 
# =============================================================================
# Extension modules
# =============================================================================
from optvl import OVLSolver
from optvl.utils.ffd_utils import write_FFD_file

try:
    from pygeo import DVGeometry
    HAS_PYGEO = True
except ImportError:
    HAS_PYGEO = False

# Add geom_files to path for importing
base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, "..", "geom_files")
sys.path.insert(0, geom_dir)

from wing_mesh import mesh  # (nx, ny, 3) numpy array

# ---------------------------------------------------------------------------
# Path setup — reuse the same mesh used in test_mesh_input.py
# ---------------------------------------------------------------------------
base_dir = os.path.dirname(os.path.abspath(__file__))
geom_dir = os.path.join(base_dir, '..', 'geom_files')
sys.path.insert(0, geom_dir)
 

surf = {
    "Wing": {
        # General
        "component": np.int32(1),  # logical surface component index (for grouping interacting surfaces, see AVL manual)
        "yduplicate": np.float64(0.0),  # surface is duplicated over the ysymm plane
        "claf": 1.0,  # CL alpha (dCL/da) scaling factor per section (provide a single entry and OptVL applies to all strips, otherwise provide a vector corresponding to each strip)

        # Geometry
        "scale": np.array(
            [1.0, 1.0, 1.0], dtype=np.float64
        ),  # scaling factors applied to all x,y,z coordinates (chords arealso scaled by Xscale)
        "translate": np.array(
            [0.0, 0.0, 0.0], dtype=np.float64
        ),  # offset added on to all X,Y,Z values in this surface
        "angle": np.float64(0.0),  # offset added on to the Ainc values for all the defining sections in this surface
        "aincs": np.ones(mesh.shape[1]), # incidence angle vector (provide a single entry and OptVL applies to all strips, otherwise provide a vector corresponding to each strip)

        # Geometry: Mesh
        "mesh": np.float64(mesh), # (nx,ny,3) numpy array containing mesh coordinates
        "flatten mesh": True, # True by default so can be turned off or just excluded (not recommended)
        
        # Control Surface Specification
        "control_assignments": {
            "flap" : {"assignment":np.arange(0,mesh.shape[1]),
                      "xhinged": 0.8, # x/c location of hinge
                      "vhinged": np.zeros(3), # vector giving hinge axis about which surface rotates
                      "gaind": 1.0, # control surface gain
                      "refld": 1.0  # control surface reflection, sign of deflection for duplicated surface
                      }
        },

        # Design Variables (AVL) Specification
        "design_var_assignments": {
            "des" : {"assignment":np.arange(0,mesh.shape[1]),
                     "gaing":1.0}
        },
    }
}

geom = {
    "title": "Aircraft",
    "mach": np.float64(0.0),
    "iysym": np.int32(0),
    "izsym": np.int32(0),
    "zsym": np.float64(0.0),
    "Sref": np.float64(10.0),
    "Cref": np.float64(1.0),
    "Bref": np.float64(10.0),
    "XYZref": np.array([0.25, 0, 0],dtype=np.float64),
    "CDp": np.float64(0.0),
    "surfaces": surf,
    # Global Control and DV info
    "dname": ["flap"],  # Name of control input for each corresonding index
    "gname": ["des"],  # Name of design var for each corresonding index
}

@unittest.skipUnless(HAS_PYGEO, "pygeo is not installed — skipping FFD tests")
class TestFFDIntegration(unittest.TestCase):
    """Tests for the pygeo FFD integration in OVLSolver
 
    Each test creates its own OVLSolver instance so that solver state is
    isolated between test cases.
    """
 
    # ------------------------------------------------------------------
    # Shared FFD file (written once for the whole test class)
    # ------------------------------------------------------------------
    ffd_path = os.path.join(geom_dir, "wing_test_ffd.xyz")
 
    def setUp(self):
        self._make_solver()
        self._make_dvgeo()

    def tearDown(self):
        # Get the memory usage of the current process using psutil
        process = psutil.Process()
        mb_memory = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
        print(f"{self.id():80} Memory usage: {mb_memory:.2f} MB")

    def _make_solver(self):
        """Sets a freshly initialised OVLSolver with the mesh geometry."""
        self.ovl_solver =  OVLSolver(input_dict=copy.deepcopy(geom))
 
    def _make_dvgeo(self):
        """Sets DVGeometry object built from the class-level FFD file to OptVL"""
        self.DVGeo = DVGeometry(self.ffd_path)
        self.ovl_solver.set_DVGeo(self.DVGeo)
        # self.ovl_solver.update_DVGeo()

    def finite_dif(self, dvgeo_seeds, step=1e-7):

        # Loop over all surfaces
        for surface in self.ovl_solver.unique_surface_names:
            # Get the pointset name
            idx_surf = self.ovl_solver.get_surface_index(surf_name=surface)
            if idx_surf in self.ovl_solver.point_sets.keys():
                point_set_name = self.ovl_solver.point_sets[idx_surf]
            else:
                continue # This surface doesn't have a pointset, skip it

            # Apply the FD step to the DV seeds
            for dv in dvgeo_seeds[surface].keys():
                # current design var values
                currentDV = self.ovl_solver.DVGeo.getValues()[dv]

                # Set the updated DVs
                self.ovl_solver.DVGeo.setDesignVars({dv: currentDV + dvgeo_seeds[surface][dv]*step})

            # get mesh size
            nx = self.ovl_solver.avl.SURF_GEOM_I.NVC[idx_surf] + 1
            ny = self.ovl_solver.avl.SURF_GEOM_I.NVS[idx_surf] + 1

            # Compute the perturbed mesh values
            coords = self.ovl_solver.DVGeo.update(point_set_name)
            mesh_pertub = copy.deepcopy(coords.reshape((ny,nx,3)).transpose((1,0,2)))

            self.ovl_solver.set_mesh(idx_surf,mesh_pertub)

        self.ovl_solver.avl.update_surfaces()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.exec_rhs()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.velsum()
        self.ovl_solver.avl.aero()
        # self.ovl_solver.execute_run()
        coef_data_peturb = self.ovl_solver.get_total_forces()
        consurf_derivs_peturb = self.ovl_solver.get_control_stab_derivs()
        stab_deriv_derivs_peturb = self.ovl_solver.get_stab_derivs()
        body_axis_deriv_petrub = self.ovl_solver.get_body_axis_derivs()
        body_forces_peturb = self.ovl_solver.get_body_forces()

        for surface in self.ovl_solver.unique_surface_names:
            # Get the pointset name
            idx_surf = self.ovl_solver.get_surface_index(surf_name=surface)
            if idx_surf in self.ovl_solver.point_sets.keys():
                point_set_name = self.ovl_solver.point_sets[idx_surf]
            else:
                continue # This surface doesn't have a pointset, skip it

            # Apply the FD step to the DV seeds
            for dv in dvgeo_seeds[surface].keys():
                # current design var values
                currentDV = self.ovl_solver.DVGeo.getValues()[dv]

                # Set the updated DVs
                self.ovl_solver.DVGeo.setDesignVars({dv: currentDV - dvgeo_seeds[surface][dv]*step})

            # get mesh size
            nx = self.ovl_solver.avl.SURF_GEOM_I.NVC[idx_surf] + 1
            ny = self.ovl_solver.avl.SURF_GEOM_I.NVS[idx_surf] + 1

            # Compute the perturbed mesh values
            coords = self.ovl_solver.DVGeo.update(point_set_name)
            mesh_orig = copy.deepcopy(coords.reshape((ny,nx,3)).transpose((1,0,2)))

            self.ovl_solver.set_mesh(idx_surf,mesh_orig)

        self.ovl_solver.avl.update_surfaces()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.exec_rhs()
        self.ovl_solver.avl.get_res()
        self.ovl_solver.avl.velsum()
        self.ovl_solver.avl.aero()
        # self.ovl_solver.execute_run()

        coef_data = self.ovl_solver.get_total_forces()
        consurf_derivs = self.ovl_solver.get_control_stab_derivs()
        stab_deriv_derivs = self.ovl_solver.get_stab_derivs()
        body_axis_deriv = self.ovl_solver.get_body_axis_derivs()
        body_forces = self.ovl_solver.get_body_forces()
        
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

    def test_mesh_deformation_propagates(self):
        self._make_solver()
        self._make_dvgeo()

        idx_surf = self.ovl_solver.get_surface_index("Wing")
        mesh_orig = copy.deepcopy(self.ovl_solver.get_mesh(idx_surf))  # (nx, ny, 3)
        DVGeo = self.ovl_solver.DVGeo
        
        DVGeo.addLocalDV("local_z", lower=-1.0, upper=1.0, axis="z", scale=1.0)
        dz = 0.1
 
        # Perturb all local z DVs by dz
        current_dvs = DVGeo.getValues()
        for key in current_dvs:
            current_dvs[key] = current_dvs[key] + dz
        DVGeo.setDesignVars(current_dvs)
 
        self.ovl_solver.update_DVGeo()
        mesh_perturb = self.ovl_solver.get_mesh(idx_surf)  # (nx, ny, 3)
 
        # The z-coordinates should have shifted by approximately dz
        dz_actual = mesh_perturb[:, :, 2] - mesh_orig[:, :, 2]
        np.testing.assert_allclose(
            dz_actual, dz,
            atol=1e-3,
            err_msg="Z-coordinates of the mesh should shift by the applied FFD dz offset",
        )
 
        # Y-coordinates should be unchanged (we only moved z)
        dy_actual = mesh_perturb[:, :, 1] - mesh_orig[:, :, 1]
        np.testing.assert_allclose(
            dy_actual, 0.0,
            atol=1e-10,
            err_msg="Y-coordinates should not change when only a Z FFD deformation is applied",
        )


    def test_totals(self):

        self._make_solver()
        self._make_dvgeo()

        # Add some DVs
        nrefaxpts = self.ovl_solver.DVGeo.addRefAxis("c4", xFraction=0.25, alignIndex="k")

        def twist(val, geo):
            for i in range(nrefaxpts):
                geo.rot_y["c4"].coef[i] = val[i]

        def sweep(val, geo):
            # the extractCoef method gets the unperturbed ref axis control points
            C = geo.extractCoef("c4")
            C_orig = C.copy()
            # we will sweep the wing about the first point in the ref axis
            sweep_ref_pt = C_orig[0, :]

            theta = -val[0] * np.pi / 180
            # rot_mtx = np.array([[np.cos(theta), 0.0, -np.sin(theta)], [0.0, 1.0, 0.0], [np.sin(theta), 0.0, np.cos(theta)]])
            rot_mtx = np.array([[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])

            # modify the control points of the ref axis
            # by applying a rotation about the first point in the x-z plane
            for i in range(nrefaxpts):
                # get the vector from each ref axis point to the wing root
                vec = C[i, :] - sweep_ref_pt
                # need to now rotate this by the sweep angle and add back the wing root loc
                C[i, :] = sweep_ref_pt + rot_mtx @ vec
            # use the restoreCoef method to put the control points back in the right place
            geo.restoreCoef(C, "c4")

        self.ovl_solver.DVGeo.addGlobalDV("twist", func=twist, value=np.zeros(nrefaxpts), lower=-10, upper=10, scale=0.05)
        self.ovl_solver.DVGeo.addGlobalDV("sweep", func=sweep, value=0.0, lower=0, upper=45, scale=0.05)
        self.ovl_solver.DVGeo.addLocalDV("shape", lower=-0.25, upper=0.25, axis="z", scale=1.0)

        # Execute the solver
        self.ovl_solver.update_DVGeo()
        self.ovl_solver.set_variable("alpha", 5.0)
        self.ovl_solver.set_variable("beta", 0.0)
        self.ovl_solver.set_parameter("Mach", 0.0)
        self.ovl_solver.execute_run()

        # compare the analytical gradients with finite difference for dvgeo
        surf_key = list(self.ovl_solver.surf_mesh_to_fort_var.keys())[0]
        dvgeo_vars = self.ovl_solver.DVGeo.getVarNames()
        cs_names = self.ovl_solver.get_control_names()

        consurf_vars = []
        for func_key in self.ovl_solver.case_derivs_to_fort_var:
            consurf_vars.append(self.ovl_solver._get_deriv_key(cs_names[0], func_key))

        func_vars = self.ovl_solver.case_var_to_fort_var
        stab_derivs = self.ovl_solver.case_stab_derivs_to_fort_var
        body_axis_derivs = self.ovl_solver.case_body_derivs_to_fort_var

        sens = self.ovl_solver.execute_run_sensitivities(
            func_vars,
            consurf_derivs=consurf_vars,
            stab_derivs=stab_derivs,
            body_axis_derivs=body_axis_derivs,
            print_timings=False,
        )

        # for con_key in self.ovl_solver.con_var_to_fort_var:
        sens_FD = {}
        for surf_key in self.ovl_solver.surf_geom_to_fort_var:
            sens_FD[surf_key] = {}
            for dvgeo_var_key in dvgeo_vars:
                # arr = self.ovl_solver.get_mesh(self.ovl_solver.get_surface_index(surf_key)).reshape(-1,3)
                arr = self.ovl_solver.DVGeo.getValues()[dvgeo_var_key]
                np.random.seed(arr.size)
                rand_arr = np.random.rand(*arr.shape)
                rand_arr /= np.linalg.norm(rand_arr)

                func_seeds, consurf_deriv_seeds, stab_derivs_seeds, body_axis_derivs_seeds = self.finite_dif({surf_key: {dvgeo_var_key: rand_arr}}, step=1.0e-7)

                for func_key in func_vars:
                    dvgeo_var_dot = np.sum(sens[func_key][surf_key][dvgeo_var_key] * rand_arr)
                    func_dot = func_seeds[func_key]

                    rel_err = np.abs(dvgeo_var_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key:5} wrt {surf_key}:{dvgeo_var_key:10} | AD:{dvgeo_var_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )
                    tol = 1e-7
                    if np.abs(dvgeo_var_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            atol=1e-4,
                            err_msg=f"{func_key:5} wrt {surf_key}:{dvgeo_var_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            rtol=5e-3,
                            err_msg=f"{func_key:5} wrt {surf_key}:{dvgeo_var_key:10}",
                        )

                for func_key in consurf_vars:
                    # for cs_key in consurf_vars[func_key]:
                    dvgeo_var_dot = np.sum(sens[func_key][surf_key][dvgeo_var_key] * rand_arr)
                    func_dot = consurf_deriv_seeds[func_key]

                    # rel_err = np.abs(dvgeo_var_dot - func_dot) / np.abs(func_dot + 1e-20)
                    # print(
                    #     f"{func_key} wrt {surf_key}:{mesh_key:10} | AD:{dvgeo_var_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )

                    tol = 1e-8
                    if np.abs(dvgeo_var_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            atol=1e-4,
                            err_msg=f"{func_key} wrt {surf_key}:{dvgeo_var_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            rtol=6e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{dvgeo_var_key:10}",
                        )

                for func_key in stab_derivs_seeds:
                    if func_key == "spiral parameter" or func_key == "lateral parameter":
                        continue # These derivatives blow up in different way for both adjoint and FD
                    dvgeo_var_dot = np.sum(sens[func_key][surf_key][dvgeo_var_key] * rand_arr)
                    func_dot = stab_derivs_seeds[func_key]

                    rel_err = np.abs(dvgeo_var_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key}  wrt {surf_key}:{dvgeo_var_key:10} | AD:{dvgeo_var_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )

                    tol = 5e-7
                    if np.abs(dvgeo_var_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            atol=5e-9,
                            err_msg=f"{func_key} wrt {surf_key}:{dvgeo_var_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            rtol=6e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{dvgeo_var_key:10}",
                        )

                for func_key in body_axis_derivs_seeds:
                    dvgeo_var_dot = np.sum(sens[func_key][surf_key][dvgeo_var_key] * rand_arr)
                    func_dot = body_axis_derivs_seeds[func_key]

                    rel_err = np.abs(dvgeo_var_dot - func_dot) / np.abs(func_dot + 1e-20)

                    # print(
                    #     f"{func_key}  wrt {surf_key}:{dvgeo_var_key:10} | AD:{dvgeo_var_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    # )

                    tol = 1e-6
                    if np.abs(dvgeo_var_dot) < tol or np.abs(func_dot) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            atol=5e-8,
                            err_msg=f"{func_key} wrt {surf_key}:{dvgeo_var_key:10}",
                        )
                    else:
                        np.testing.assert_allclose(
                            dvgeo_var_dot,
                            func_dot,
                            rtol=6e-3,
                            err_msg=f"{func_key} wrt {surf_key}:{dvgeo_var_key:10}",
                        )


if __name__ == "__main__":
    unittest.main()