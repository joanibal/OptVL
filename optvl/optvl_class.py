"""
OptVL

OptVL is a wrapper for Mark Drela's AVL code. The purpose of this
class is to provide an easy to use wrapper for avl for intergration
into other projects.

Developers:
-----------
- Josh Anibal (JLA)

History
-------
    v. 0.0 - Initial Class Creation (JLA, 08 2016)
    v. 0.1 - Added CDFF and sectional properties for all surfaces (JLA, 12 2017)
    v. 1.0 - major refactor (I learned a lot more programming in 7 years) (JLA, 02 2023)
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import copy
from pprint import pprint
from typing import Dict, List, Tuple, Union, Any
import warnings
import glob
from typing import Optional
import platform

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np


# =============================================================================
# Extension modules
# =============================================================================
from . import MExt


class AVLSolver(object):

    # these at technically parameters, but they are also specified as contraints
    # These are not included in the derivatives but you can set and get them still
    # In the code these are only used to save the state for mode analysis
    # regardless they can be useful for getting the values after a trim solve
    state_param_idx_dict = {
        "alpha": 0, 
        "beta": 1,
        "roll rate": 2,
        "pitch rate": 3,
        "yaw rate": 4,
        # Not a supported way to set CL as a contraint or get it as an output
        # "CL": 5, 
    }

    param_idx_dict = {
        # these aero state parameters
        "CD0": 6,
        "bank": 7,
        "elevation": 8,
        "heading": 9,
        "Mach": 10,
        "velocity": 11,
        "density": 12,
        "grav.acc.": 13,
        "turn rad.": 14,
        "load fac.": 15,
        "X cg": 16,
        "Y cg": 17,
        "Z cg": 18,
        "mass": 19,
        "Ixx": 20,
        "Iyy": 21,
        "Izz": 22,
        "Ixy": 23,
        "Iyz": 24,
        "Izx": 25,
        "visc CL_a": 26,
        "visc CL_u": 27,
        "visc CM_a": 28,
        "visc CM_u": 29,
    }
    
    # fmt: off
    # This dict has the following structure:
    # python key: [common block name, fortran varaiable name]
    case_var_to_fort_var = {
        # lift and drag from surface integration (wind frame)
        "CL": ["CASE_R", "CLTOT"],
        "CD": ["CASE_R", "CDTOT"],
        "CDv": ["CASE_R", "CDVTOT"], # viscous drag
        
        # lift and drag calculated from farfield integration
        "CLff": ["CASE_R", "CLFF"],
        "CYff": ["CASE_R", "CYFF"],
        "CDi": ["CASE_R", "CDFF"], # induced drag
        
        # non-dimensionalized forces 
        "CX": ["CASE_R", "CXTOT"], 
        "CY": ["CASE_R", "CYTOT"], 
        "CZ": ["CASE_R", "CZTOT"],   
        
        # non-dimensionalized moments (body frame)
        "CR BA": ["CASE_R", "CRTOT"],
        "CM": ["CASE_R", "CMTOT"],
        "CN BA": ["CASE_R", "CNTOT"],

        # non-dimensionalized moments (stablity frame)
        "CR SA": ["CASE_R", "CRSAX"],
        # "CM SA": ["CASE_R", "CMSAX"], # This is the same in both frames
        "CN SA": ["CASE_R", "CNSAX"],
        
        # spanwise efficiency
        "e": ["CASE_R", "SPANEF"],
    }
    ref_var_to_fort_var = {
        "Sref": ["CASE_R", "SREF"],
        "Cref": ["CASE_R", "CREF"],
        "Bref": ["CASE_R", "BREF"],
    }
    
    case_derivs_to_fort_var = {
        # derivative of coefficents wrt control surface deflections
        "CL": ["CASE_R", "CLTOT_D"],
        "CD": ["CASE_R", "CDTOT_D"],
        "CX": ["CASE_R", "CXTOT_D"],
        "CY": ["CASE_R", "CYTOT_D"],
        "CZ": ["CASE_R", "CZTOT_D"],
        "CR": ["CASE_R", "CRTOT_D"],
        "CM": ["CASE_R", "CMTOT_D"],
        "CN": ["CASE_R", "CNTOT_D"],
    }
  
    # This dict has the following structure:
    # python key: [common block name, fortran varaiable name]
    case_surf_var_to_fort_var = {
        # surface contributions to total lift and drag from surface integration (wind frame)
        "CL": ["SURF_R", "CLSURF"],
        "CD": ["SURF_R", "CDSURF"],
        "CDv": ["SURF_R", "CDVSURF"],  # viscous drag
        
        # non-dimensionalized forces 
        "CX": ["SURF_R", "CXSURF"], 
        "CY": ["SURF_R", "CYSURF"], 
        "CZ": ["SURF_R", "CZSURF"],   
        
        # non-dimensionalized moments (body frame)
        "CR": ["SURF_R", "CRSURF"],
        "CM": ["SURF_R", "CMSURF"],
        "CN": ["SURF_R", "CNSURF"],
        
        # forces non-dimentionalized by surface quantities
        # uses surface area instead of sref and takes moments about leading edge
        "CL surf" : ["SURF_R", "CL_SRF"],
        "CD surf" : ["SURF_R", "CD_SRF"],
        "CMLE surf" : ["SURF_R", "CMLE_SRF"],
        
        #TODO: add CF_SRF(3,NFMAX), CM_SRF(3,NFMAX)
    }
    

    body_geom_to_fort_var = {
        "scale": ["BODY_GEOM_R", "XYZSCAL_B"],
        "translate": ["BODY_GEOM_R", "XYZTRAN_B"],
        "yduplicate": ["BODY_GEOM_R", "YDUPL_B"],
        "bfile": ["CASE_C", "BFILES"],
        "nvb": ["BODY_GEOM_I", "NVB"],
        "bspace": ["BODY_GEOM_R", "BSPACE"],

    }
    
    # fmt: on

    ad_suffix = "_DIFF"
    
    NUMAX = 6

    def __init__(self, geo_file=None, mass_file=None, debug=False, timing=False):

        if timing:
            start_time = time.time()

        # MExt is important for creating multiple instances of the AVL solver that do not share memory
        # It is very gross, but I cannot figure out a better way (maybe use install_name_tool to change the dynamic library path to absolute).
        # increment this counter for the hours you wasted on trying find a better way
        # 7 hours

        module_dir = os.path.dirname(os.path.realpath(__file__))
        module_name = os.path.basename(module_dir)
        
        if platform.system() == "Windows":
            #HACK
            avl_lib_so_file = glob.glob(os.path.join(module_dir, "libavl*.pyd"))[0]
        else:
            avl_lib_so_file = glob.glob(os.path.join(module_dir, "libavl*.so"))[0]
            
        # # get just the file name
        avl_lib_so_file = os.path.basename(avl_lib_so_file)
        self.avl = MExt.MExt("libavl", module_name, "optvl", lib_so_file=avl_lib_so_file, debug=debug)._module

        # this way doesn't work with mulitple isntances fo AVLSolver
        # from . import libavl
        # self.avl = libavl

        if not (geo_file is None):
            try:
                # check to make sure files exist
                file = geo_file
                f = open(geo_file, "r")
                f.close()

                if not (mass_file is None):
                    file = mass_file
                    f = open(mass_file, "r")
                    f.close()
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not open the file '{file}' from python. This is usually an issue with the specified file path"
                )

            self.avl.avl()
            if debug:
                self.set_avl_fort_arr("CASE_L", "LVERBOSE", True)

            if timing:
                self.set_avl_fort_arr("CASE_L", "LTIMING", True)

            self.avl.loadgeo(geo_file)

            if mass_file is not None:
                self.avl.loadmass(mass_file)

        else:
            raise ValueError("neither a geometry file or aircraft object was given")
        
        # todo store the default dict somewhere else
        # the control surface contraints get added to this array in the __init__
        self.conval_idx_dict = {
            "alpha": 0,
            "beta": 1,
            "roll rate": 2,
            "pitch rate": 3,
            "yaw rate": 4,
            "CL": 5,
            "CY": 6,
            "CR BA": 7,
            "CM": 8,
            "CR": 9,
        }
        
        # control surfaces added in __init__
        #TODO: the keys of this dict aren't used
        self.con_var_to_fort_var = {
            "alpha": ["CASE_R", "ALFA"],
            "beta": ["CASE_R", "BETA"],
        }



        control_names = self.get_control_names()
        self.control_variables = {}
        for idx_c_var, c_name in enumerate(control_names):
            self.control_variables[c_name] = f"D{idx_c_var+1}"
            
        # set control surface constraint indecies in to con val dict
        idx_control_start = np.max([x for x in self.conval_idx_dict.values()]) + 1
        for idx_c_var, c_name in enumerate(control_names):
            self.conval_idx_dict[c_name] = idx_control_start + idx_c_var
            self.con_var_to_fort_var[c_name] = ["CASE_R", "DELCON"]
            
        var_to_suffix = {
            "alpha": "AL",
            "beta": "BE",
            "roll rate": "RX",
            "pitch rate": "RY",
            "yaw rate": "RZ",
            }
        self.case_stab_derivs_to_fort_var = {}
        for func in ["CL", "CD", "CY", "CR", "CM", "CN"]:
            for var in var_to_suffix:
                self.case_stab_derivs_to_fort_var[self.get_deriv_key(var, func)] = ["CASE_R", f"{func}TOT_{var_to_suffix[var]}"]
                
        #  the case parameters are stored in a 1d array,
        # these indices correspond to the position of each parameter in that arra
        self._init_surf_data()

        # set the default solver tolerance
        self.set_avl_fort_arr('CASE_R', 'EXEC_TOL', 2e-5)

        if timing:
            print(f"AVL init took {time.time() - start_time} seconds")

    def _init_surf_data(self):
        self.surf_geom_to_fort_var = {}
        surf_names = self.get_surface_names()

        self.surf_pannel_to_fort_var = {}
        self.con_surf_to_fort_var = {}

        for idx_surf, surf_name in enumerate(surf_names):
            idx_surf = surf_names.index(surf_name)
            slice_idx_surf = (idx_surf,)
            slice_surf_all = (idx_surf, slice(None))

            # only set unduplicated sufaces
            if self.get_avl_fort_arr("SURF_I", "IMAGS", slicer=slice_idx_surf) < 0:
                # this is a duplicated surface, skip it
                continue

            num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC", slicer=slice_idx_surf)

            slice_surf_secs = (idx_surf, slice(None, num_sec))
            slice_surf_secs_all = slice_surf_secs + (slice(None),)

            nasec = self.get_avl_fort_arr("SURF_GEOM_I", "NASEC", slicer=slice_surf_secs)
            # I don't see a case in the code where nasec is not the same for all sections
            # but just in case, we'll do a test and throw an error if not
            if np.unique(nasec).size != 1:
                raise RuntimeError("nasec is not the same for all sections")

            nasec = nasec[0]

            slice_surf_secs_nasec = slice_surf_secs + (slice(None, nasec),)

            self.surf_geom_to_fort_var[surf_name] = {
                "scale": ["SURF_GEOM_R", "XYZSCAL", slice_surf_all],
                "translate": ["SURF_GEOM_R", "XYZTRAN", slice_surf_all],
                "angle": ["SURF_GEOM_R", "ADDINC", slice_idx_surf],
                "xyzles": ["SURF_GEOM_R", "XYZLES", slice_surf_secs_all],
                "chords": ["SURF_GEOM_R", "CHORDS", slice_surf_secs],
                "aincs": ["SURF_GEOM_R", "AINCS", slice_surf_secs],
                "xasec": ["SURF_GEOM_R", "XASEC", slice_surf_secs_nasec],
                "sasec": ["SURF_GEOM_R", "SASEC", slice_surf_secs_nasec],
                "tasec": ["SURF_GEOM_R", "TASEC", slice_surf_secs_nasec],
                "clcdsec": ["SURF_GEOM_R", "CLCDSEC", slice_surf_secs_all],
                "claf": ["SURF_GEOM_R", "CLAF", slice_surf_secs],
            }
            self.surf_pannel_to_fort_var[surf_name] = {
                "nchordwise": ["SURF_GEOM_I", "NVC", slice_idx_surf],
                "cspace": ["SURF_GEOM_R", "CSPACE", slice_idx_surf],
                "nspan": ["SURF_GEOM_I", "NVS", slice_idx_surf],
                "sspace": ["SURF_GEOM_R", "SSPACE", slice_idx_surf],
                "sspaces": ["SURF_GEOM_R", "SSPACES", slice_surf_secs],
                "nspans": ["SURF_GEOM_I", "NSPANS", slice_surf_secs],
                "yduplicate": ["SURF_GEOM_R", "YDUPL", slice_idx_surf],
                "use surface spacing": ["SURF_GEOM_L", "LSURFSPACING", slice_idx_surf],
                "component": ["SURF_I", "LSCOMP", slice_idx_surf],
            }

            icontd_slices = []
            idestd_slices = []
            xhinged_slices = []
            vhinged_slices = []
            gaind_slices = []
            refld_slices = []
            gaing_slices = []

            for idx_sec in range(num_sec):
                slice_surf = (idx_surf, idx_sec)

                num_con_surf = self.get_avl_fort_arr("SURF_GEOM_I", "NSCON", slicer=slice_surf)
                slice_surf_con = slice_surf + (slice(None, num_con_surf),)
                slice_surf_con_hinge = slice_surf_con + (slice(None),)

                num_des_var = self.get_avl_fort_arr("SURF_GEOM_I", "NSDES", slicer=slice_surf)
                slice_surf_des_var = slice_surf + (slice(None, num_des_var),)

                icontd_slices.append(slice_surf_con)
                xhinged_slices.append(slice_surf_con)
                vhinged_slices.append(slice_surf_con_hinge)
                gaind_slices.append(slice_surf_con)
                refld_slices.append(slice_surf_con)
                idestd_slices.append(slice_surf_des_var)
                gaing_slices.append(slice_surf_des_var)

            self.con_surf_to_fort_var[surf_name] = {
                "icontd": ["SURF_GEOM_I", "ICONTD", icontd_slices],
                "xhinged": ["SURF_GEOM_R", "XHINGED", xhinged_slices],
                "vhinged": ["SURF_GEOM_R", "VHINGED", vhinged_slices],
                "gaind": ["SURF_GEOM_R", "GAIND", gaind_slices],
                "refld": ["SURF_GEOM_R", "REFLD", refld_slices],
                "idestd": ["SURF_GEOM_I", "IDESTD", idestd_slices],
                "gaing": ["SURF_GEOM_R", "GAING", gaing_slices],
            }

# region -- analysis api
    def execute_run(self, tol=0.00002):
        # run the analysis (equivalent to the avl command `x` in the oper menu)
        self.set_avl_fort_arr('CASE_R', 'EXEC_TOL', tol)
        self.avl.oper()

    def set_constraint(self, var, val, con_var=None):
        avl_variables = {
            "alpha": "A",
            "beta": "B",
            "roll rate": "R",
            "pitch rate": "P",
            "yaw rate": "Y",
        }

        avl_con_variables = copy.deepcopy(avl_variables)
        avl_con_variables.update(
            {
                "CL": "C",
                "CY": "S",
                "Cl roll moment": "RM",
                "Cm pitch moment": "PM",
                "Cn yaw moment": "YM",
            }
        )

        if var in avl_variables:
            # save the name of the avl_var
            avl_var = avl_variables[var]
        elif var in self.control_variables.keys():
            avl_var = self.control_variables[var]
        elif var in self.control_variables.values():
            avl_var = var
        else:
            raise ValueError(
                f"specified variable `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}. Constraints that must be implicitly satisfied (such as `CL`) are set with `add_trim_constraint`."
            )

        if con_var is None:
            avl_con_var = avl_var
        elif con_var in avl_con_variables:
            avl_con_var = avl_con_variables[con_var]
        elif con_var in self.control_variables.keys():
            avl_con_var = self.control_variables[con_var]
        elif con_var in self.control_variables.values():
            avl_con_var = con_var
        else:
            raise ValueError(
                f"specified contraint variable `{con_var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}."
            )

        # check that the type of val is correct
        if not isinstance(val, (int, float, np.floating, np.integer)):
            raise TypeError(f"contraint value must be a int or float for contraint {var}. Got {type(val)}")

        self.avl.conset(avl_var, f"{avl_con_var} {val} \n")

    def set_trim_condition(self, variable, val):
        options = {
            "bankAng": ["B"],
            "CL": ["C"],
            "velocity": ["V"],
            "mass": ["M"],
            "dens": ["D"],
            "G": ["G"],
            "X cg": ["X"],
            "Y cg": ["Y"],
            "Z cg": ["Z"],
        }

        if not (variable in options):
            raise ValueError(
                f"constraint variable `{variable}` not a valid option. Must be one of the following {[key for key in options]} "
            )

        self.avl.trmset("C1", "1 ", options[variable][0], (str(val) + "  \n"))

    def get_total_forces(self) -> Dict[str, float]:
        """Get the aerodynamic data for the last run case and return it as a dictionary.

        Returns:
            Dict[str, float]: Dictionary of aerodynamic data. The keys the aerodyanmic coefficients.
        """

        total_data = {}

        for key, avl_key in self.case_var_to_fort_var.items():
            val = self.get_avl_fort_arr(*avl_key)
            # [()] because all the data is stored as a ndarray.
            # for scalars this results in a 0-d array.
            # It is easier to work with floats so we extract the value with [()]
            total_data[key] = val[()]

        return total_data

    def get_control_stab_derivs(self) -> Dict[str, Dict[str, float]]:
        deriv_data = {}

        control_names = self.get_control_names()

        for key, avl_key in self.case_derivs_to_fort_var.items():
            slicer = slice(0, len(control_names))
            val_arr = self.get_avl_fort_arr(*avl_key, slicer=slicer)
            for idx_control, val in enumerate(val_arr):
                control = control_names[idx_control]
                deriv_data[self.get_deriv_key(control, key)] = val[()]

        return deriv_data

    def get_stab_derivs(self) -> Dict[str, Dict[str, float]]:
        deriv_data = {}

        for func_key, avl_key in self.case_stab_derivs_to_fort_var.items():

            val_arr = self.get_avl_fort_arr(*avl_key)
            deriv_data[func_key] = val_arr[()]

        return deriv_data
    
    def get_reference_data(self) -> Dict[str, float]:
        ref_data = {}
        
        for key, avl_key in self.ref_var_to_fort_var.items():
            ref_data[key] = self.get_avl_fort_arr(*avl_key)[()]
            
        return ref_data
        
    def set_reference_data(self, ref_data: Dict[str, float]) -> None:
        for key, val in ref_data.items():
            avl_key = self.ref_var_to_fort_var[key]
            self.set_avl_fort_arr(*avl_key, val)
            
        return ref_data

    def get_avl_fort_arr(self, common_block, variable, slicer=None):
        # this had to be split up into two steps to work

        # get the corresponding common block object.
        # it must be lowercase because of f2py
        common_block = getattr(self.avl, common_block.upper())

        # get the value of the variable from the common block
        val = getattr(common_block, variable.upper())

        # convert from fortran ordering to c ordering
        val = val.ravel(order="F").reshape(val.shape[::-1], order="C")

        # Apply slicer if provided
        if slicer is not None:
            val = val[slicer]

        return val

    def set_avl_fort_arr(self, common_block, variable, val, slicer=None) -> None:
        # convert from fortran ordering to c ordering
        if isinstance(val, np.ndarray):
            val = val.ravel(order="C").reshape(val.shape[::-1], order="F")

        # this had to be split up into two steps to work
        # get the corresponding common block object.
        # it must be lowercase because of f2py
        common_block_obj = getattr(self.avl, common_block.upper())

        # get the value of the variable from the common block
        if slicer is None:
            setattr(common_block_obj, variable.upper(), val)
        else:
            # flip the order of the slicer to match the cordinates of the val
            new_slicer = slicer[::-1]

            original_val = getattr(common_block_obj, variable.upper())
            original_val[new_slicer] = val
            setattr(common_block_obj, variable.upper(), original_val)

        return

    def get_surface_forces(self) -> Dict[str, Dict[str, float]]:
        surf_names = self.get_surface_names()

        # add a dictionary for each surface that will be filled later
        surf_data = {}
        for surf in surf_names:
            surf_data[surf] = {}

        for key, avl_key in self.case_surf_var_to_fort_var.items():
            vals = self.get_avl_fort_arr(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(surf_names):
                surf_data[surf_name][key] = vals[idx_surf]

        return surf_data

    def get_parameter(self, param_key: str) -> float:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")

        # the key could be in one of two dicts
        if param_key in self.param_idx_dict:
            idx_param = self.param_idx_dict[param_key]
        elif param_key in self.state_param_idx_dict:
            idx_param = self.state_param_idx_dict[param_key]
        else:
            raise ValueError(f"param '{param_key}' not in possilbe list\n" 
                             f"{[k for k in self.param_idx_dict] + [k for k in self.state_param_idx_dict]}")
        # [0] because optvl only supports 1 run case
        param_val = parvals[0][idx_param]

        return param_val

    def get_constraint(self, con_key: str) -> float:
        """ """
        convals = self.get_avl_fort_arr("CASE_R", "CONVAL")

        # [0] because optvl only supports 1 run case
        con_val = convals[0][self.conval_idx_dict[con_key]]

        return con_val

    def set_parameter(self, param_key: str, param_val: float) -> None:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        # warn the user that alpha, beta,
        if param_key is ["alpha", "beta", "pb/2V", "qc/2V", "rb/2V", "CL"]:
            raise ValueError(
                "alpha, beta, pb/2V, qc/2V, rb/2V, and CL are not allowed to be set,\n\
                             they are calculated during each run based on the contraints. to specify\n\
                             one of these values use the add_contraint method."
            )

        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")
        # [0] because optvl only supports 1 run case
        parvals[0][self.param_idx_dict[param_key]] = param_val

        self.set_avl_fort_arr("CASE_R", "PARVAL", parvals)

        # (1) here because we want to set the first runcase with fortran indexing (the only one)
        self.avl.set_params(1)

    def get_control_deflections(self) -> Dict[str, float]:
        control_surfaces = self.get_control_names()

        def_arr = copy.deepcopy(self.get_avl_fort_arr("CASE_R", "DELCON"))

        def_dict = {}
        for idx_con, con_surf in enumerate(control_surfaces):
            def_dict[con_surf] = def_arr[idx_con]

        return def_dict

    def get_hinge_moments(self) -> Dict[str, float]:
        """
        get the hinge moments from the fortran layer and return them as a dictionary
        """
        hinge_moments = {}

        control_surfaces = self.get_control_names()
        mom_array = self.get_avl_fort_arr("CASE_R", "CHINGE")

        for idx_con, con_surf in enumerate(control_surfaces):
            hinge_moments[con_surf] = mom_array[idx_con]

        return hinge_moments

    def get_strip_forces(self) -> Dict[str, Dict[str, np.ndarray]]:
        # fmt: off
        var_to_fort_var = {
            # geometric quantities
            "chord": ["STRP_R", "CHORD"],
            "width": ["STRP_R", "WSTRIP"],
            "XYZ LE": ["STRP_R", "RLE"],  # control point leading edge coordinates
            "twist": ["STRP_R", "AINC"],
            
            # strip contributions to total lift and drag from strip integration
            "CL": ["STRP_R", "CLSTRP"],
            "CD": ["STRP_R", "CDSTRP"],
            "CDv" : ["STRP_R","CDV_LSTRP"],  # strip viscous drag in stability axes
            "downwash" : ["STRP_R","DWWAKE"],
            
            
            # strip contributions to non-dimensionalized forces 
            "CX": ["STRP_R", "CXSTRP"], 
            "CY": ["STRP_R", "CYSTRP"], 
            "CZ": ["STRP_R", "CZSTRP"],   
            
            # strip contributions to total moments (body frame)
            "CM": ["STRP_R", "CMSTRP"],
            "CN": ["STRP_R", "CNSTRP"],
            "CR": ["STRP_R", "CRSTRP"],
            
            
            # forces non-dimentionalized by strip quantities
            "CL strip" : ["STRP_R", "CL_LSTRP"],
            "CD strip" : ["STRP_R", "CD_LSTRP"],
            "CF strip" : ["STRP_R", "CF_STRP"], # forces in 3 directions
            "CM strip" : ["STRP_R", "CM_STRP"], # moments in 3 directions

            # additional forces and moments
            "CL perp" : ["STRP_R", "CLTSTRP"], # strip CL referenced to Vperp,
            "CM c/4" : ["STRP_R","CMC4"],  # strip pitching moment about c/4 and
            "CM LE" : ["STRP_R","CMLE"],  # strip pitching moment about LE vector
            "spanloading" : ["STRP_R","CNC"],   # strip spanloading 
            
        }    
        # fmt: on
        surf_names = self.get_surface_names()

        # add a dictionary for each surface that will be filled later
        strip_data = {}
        for surf in surf_names:
            strip_data[surf] = {}

        for key, avl_key in var_to_fort_var.items():
            vals = self.get_avl_fort_arr(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(surf_names):
                idx_srp_beg, idx_srp_end = self.get_surface_strip_indices(idx_surf)
                strip_data[surf_name][key] = vals[idx_srp_beg:idx_srp_end]

        cref = self.get_avl_fort_arr("CASE_R", "CREF")
        for surf_key in strip_data:
            # add sectional lift and drag
            strip_data[surf_key]["lift dist"] = strip_data[surf_key]["CL"] * strip_data[surf_key]["chord"] / cref
            strip_data[surf_key]["drag dist"] = strip_data[surf_key]["CD"] * strip_data[surf_key]["chord"] / cref

        return strip_data

    def get_surface_strip_indices(self, idx_surf):
        num_strips = np.trim_zeros(self.get_avl_fort_arr("SURF_I", "NJ"))
        idx_srp_beg = np.sum(num_strips[:idx_surf])
        idx_srp_end = np.sum(num_strips[: idx_surf + 1])

        return idx_srp_beg, idx_srp_end

# region --- modal analysis api
    def execute_eigen_mode_calc(self):
        self.avl.execute_eigenmode_calc()
    
    def get_eigenvalues(self):
        """after running an eigenmode calculation, this function will return the eigenvalues in the order used by AVL"""
        
        # get the number of "valid" eigenvalues from avl
        # [0] because optvl only supports 1 run case
        num_eigen = self.get_avl_fort_arr("CASE_I", "NEIGEN")[0]
        
        # 0 because optvl only supports 1 run case
        slicer = (0, slice(0,num_eigen))
        # get the eigenvalues from avl
        eig_vals = self.get_avl_fort_arr("CASE_Z", "EVAL", slicer=slicer)
        return eig_vals

    def get_eigenvectors(self):
        """after running an eigenmode calculation, this function will return the eigenvalues in the order used by AVL"""
        
        # get the number of "valid" eigenvalues from avl
        # [0] because optvl only supports 1 run case
        num_eigen = self.get_avl_fort_arr("CASE_I", "NEIGEN")[0]
        
        # 0 because optvl only supports 1 run case
        slicer = (0, slice(0,num_eigen), slice(None))
        eig_vecs = self.get_avl_fort_arr("CASE_Z", "EVEC", slicer=slicer)
        
        return eig_vecs
    
    def get_system_matrix(self):
        """returns the system matrix used for the eigenmode calculation"""
        
        # get the dimesion of the A matrix from the eig_vals    
        eig_vals = self.get_avl_fort_arr("CASE_Z", "EVAL")
        jemax = eig_vals.shape[1]
        asys = np.zeros((jemax,jemax), order="F")
        
        # 1 because optvl only supports 1 run case and we are using fortran base 1 indexing
        irun_case = 1
        self.avl.get_system_matrix(irun_case,asys)
        
        return asys

# region --- geometry api
    def get_control_names(self) -> List[str]:
        fort_names = self.get_avl_fort_arr("CASE_C", "DNAME")
        control_names = self._convertFortranStringArrayToList(fort_names)
        return control_names

    def get_surface_names(self, remove_dublicated=False) -> List[str]:
        """get the surface names from the geometry"""
        fort_names = self.get_avl_fort_arr("CASE_C", "STITLE")
        surf_names = self._convertFortranStringArrayToList(fort_names)

        if remove_dublicated:
            imags = self.get_avl_fort_arr("SURF_I", "IMAGS")
            unique_surf_names = []
            for idx_surf, surf_name in enumerate(surf_names):
                # get surfaces that have not been duplicated
                if imags[idx_surf] > 0:
                    unique_surf_names.append(surf_names[idx_surf])

            return unique_surf_names
        else:
            return surf_names

    def get_body_names(self, remove_dublicated=False) -> List[str]:
        """get the body names from the geometry"""
        fort_names = self.get_avl_fort_arr("CASE_C", "BTITLE")
        body_names = self._convertFortranStringArrayToList(fort_names)

        if remove_dublicated:
            # imags = self.get_avl_fort_arr("BODY_GEOM_L", "LDUPL_B")
            # print(imags)
            unique_body_names = []
            
            for body_name in body_names:
                # get bodyaces that have not been duplicated
                
                # HACK: It is best not rely on this but, this is a quick fix for 
                # bodies which I discourage people from using anyways
                if not body_name.endswith("(YDUP)"):
                    unique_body_names.append(body_name)

            return unique_body_names
        else:
            return body_names

    def get_con_surf_param(self, surf_name, idx_slice, param):
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.con_surf_to_fort_var[surf_name].keys():
            fort_var = self.con_surf_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has control surface data {self.con_surf_to_fort_var[surf_name].keys()}"
            )
        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_slice])
        return param

    def set_con_surf_param(self, surf_name, idx_slice, param, val):
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.con_surf_to_fort_var[surf_name].keys():
            fort_var = self.con_surf_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has control surface data {self.con_surf_to_fort_var[surf_name].keys()}"
            )
        # param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_slice])
        self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2][idx_slice])

    def get_surface_param(self, surf_name, param):
        # check that param is in self.surf_geom_to_fort_var
        if param in self.surf_geom_to_fort_var[surf_name].keys():
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
        elif param in self.surf_pannel_to_fort_var[surf_name].keys():
            fort_var = self.surf_pannel_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )

        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2])
        return param

    def set_surface_param(self, surf_name, param, val):
        if param in self.surf_geom_to_fort_var[surf_name].keys():
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
        elif param in self.surf_pannel_to_fort_var[surf_name].keys():
            fort_var = self.surf_pannel_to_fort_var[surf_name][param]
        elif param in self.con_surf_to_fort_var[surf_name].keys():
            # the control surface and design variables need to be handeled differently because the number at each section is variable
            pass

        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )

        self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2])

    def get_surface_params(
        self,
        include_geom: bool = True,
        include_panneling: bool = False,
        include_con_surf: bool = False,
        include_airfoils: bool = False,
    ) -> Dict[str, Dict[str, Any]]:
        """get the surface level parameters from the geometry

        geom_only: only return the geometry parameters of the surface
            - xle, yle, zle, chord, twist
        """

        surf_names = self.get_surface_names()
        unique_surf_names = self.get_surface_names(remove_dublicated=True)
        surf_data = {}

        for surf_name in unique_surf_names:
            surf_data[surf_name] = {}
            if include_geom:
                for var in self.surf_geom_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

            idx_surf = surf_names.index(surf_name)
            if include_panneling:
                # add panneling parameters if requested
                for var in self.surf_pannel_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

                if not self.get_avl_fort_arr("SURF_GEOM_L", "LDUPL")[idx_surf]:
                    surf_data[surf_name].pop("yduplicate")

            if include_con_surf:
                # add control surface parameters if requested
                for var in self.con_surf_to_fort_var[surf_name]:
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]

                    slice_data = []
                    for idx_sec in range(num_sec):
                        tmp = self.get_con_surf_param(surf_name, idx_sec, var)
                        slice_data.append(tmp)

                    surf_data[surf_name][var] = slice_data

            if include_airfoils:
                afiles = []
                num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]

                for idx_sec in range(num_sec):
                    afile = self.__decodeFortranString(self.avl.CASE_C.AFILES[idx_sec, idx_surf])
                    afiles.append(afile)

                surf_data[surf_name]["afiles"] = afiles

        return surf_data

    def set_surface_params(self, surf_data: Dict[str, Dict[str, any]]) -> None:
        """set the give surface data into the current avl object.
        ASSUMES THE CONTROL SURFACE DATA STAYS AT THE SAME LOCATION"""
        surf_names = self.get_surface_names()
        unique_surf_names = self.get_surface_names(remove_dublicated=True)

        for surf_name in surf_data:
            if surf_name not in unique_surf_names:
                raise ValueError(
                    f"""surface name, {surf_name}, not found in the current avl object."
                        Note duplicated surfaces can not be set directly.
                        Surface in file {unique_surf_names}
                        {surf_names}"""
                )

            for var in surf_data[surf_name]:
                # do not set the data this way if it is a control surface
                if var not in self.con_surf_to_fort_var[surf_name]:
                    self.set_surface_param(surf_name, var, surf_data[surf_name][var])
                else:
                    idx_surf = surf_names.index(surf_name)
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]
                    slice_data = []
                    for idx_sec in range(num_sec):
                        self.set_con_surf_param(surf_name, idx_sec, var, surf_data[surf_name][var][idx_sec])

        self.avl.update_surfaces()
    
    def get_body_params(self) ->  Dict[str, Dict[str, Any]]:
        
        body_names = self.get_body_names()
        unique_body_names = self.get_body_names(remove_dublicated=True)
        body_data = {}

        for body_name in unique_body_names:
            idx_body = body_names.index(body_name)
            body_data[body_name] = {}
            
            for var, fort_var in self.body_geom_to_fort_var.items():
                val = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=idx_body)            # 
                if var == "bfile":
                    val = self.__decodeFortranString(val)
                    
                body_data[body_name][var] = val
        
        return body_data

# region --- geometry file writing api
    def write_geom_file(self, filename: str):
        """write the current avl geometry to a file"""
        with open(filename, "w") as fid:
            # write the header
            fid.write("# generated using OptVL\n")
            self.__write_header(fid)

            surf_data = self.get_surface_params(
                include_geom=True, include_panneling=True, include_con_surf=True, include_airfoils=True
            )

            for surf_name in surf_data:
                self.__write_surface(fid, surf_name, surf_data[surf_name])
                
            body_data = self.get_body_params()
            for body_name in body_data:
                self.__write_body(fid, body_name, body_data[body_name])

    def __write_fort_vars(self, fid, common_block: str, fort_var: str, newline: bool = True) -> None:
        var = self.get_avl_fort_arr(common_block, fort_var)

        out_str = ""
        # loop over the variables list and recursively convert the variables to a string and add to the output string
        if isinstance(var, np.ndarray):
            if var.size == 1:
                out_str += str(var[()])
            else:
                out_str += " ".join([str(item) for item in var])
        else:
            out_str += str(var)
        out_str += " "

        if newline:
            out_str += "\n"

        fid.write(out_str)

    def __write_header(self, fid):
        """write the header to a file"""
        # write the name of the aircraft

        self.__write_banner(fid, "Header")
        title_array = self.get_avl_fort_arr("case_c", "title")
        title = self._convertFortranStringArrayToList(title_array)
        fid.write(f"{title}\n")

        fid.write("#Mach\n")
        self.__write_fort_vars(fid, "case_r", "mach0")

        fid.write("#IYsym   IZsym   Zsym\n")
        self.__write_fort_vars(fid, "case_i", "iysym", newline=False)
        self.__write_fort_vars(fid, "case_i", "izsym", newline=False)
        self.__write_fort_vars(fid, "case_r", "zsym")

        fid.write("#Sref    Cref    Bref\n")
        self.__write_fort_vars(fid, "case_r", "sref", newline=False)
        self.__write_fort_vars(fid, "case_r", "cref", newline=False)
        self.__write_fort_vars(fid, "case_r", "bref")

        fid.write("#Xref    Yref    Zref\n")
        self.__write_fort_vars(fid, "case_r", "XYZREF")

        fid.write("#CD0\n")
        self.__write_fort_vars(fid, "case_r", "CDREF0")

        # fid.write(f" {self.get_avl_fort_arr('case_r', 'sref')}")

    def __write_banner(self, fid, header, line_width: int = 80):
        header = " " + header + " "  # pad with spaces

        width = line_width - 1

        banner = f"#{'='*(width)}\n" f"#{header.center(width,'-')}\n" f"#{'='*(width)}\n"
        fid.write(banner)

        # ======================================================
        # ------------------- Geometry File --------------------
        # ======================================================

    def __write_body(self, fid, body_name, data):
        self.__write_banner(fid, body_name)
        fid.write(f"BODY\n")
        fid.write(f"{body_name}\n")
        fid.write(f"#N  Bspace\n")
        fid.write(f"{data['nvb']} {data['bspace']}\n")
        fid.write("SCALE\n")
        fid.write(f"{data['scale'][0]} {data['scale'][1]} {data['scale'][2]}\n")
        fid.write("TRANSLATE\n")
        fid.write(f"{data['translate'][0]} {data['translate'][1]} {data['translate'][2]}\n")
        if data['bfile'] != '':
            fid.write(f"BFILE\n")
            fid.write(f"{data['bfile']}\n")
        
    def __write_surface(self, fid, surf_name, data):
        """write a surface to a file"""
        # TODO add NACA and CLAF keyword support

        def __write_data(key_list, newline: bool = True):
            out_str = ""
            for key in key_list:
                val = data[key]

                if isinstance(val, np.ndarray):
                    if val.size == 1:
                        out_str += str(val[()])
                    else:
                        out_str += " ".join([str(item) for item in val])
                else:
                    out_str += str(val)
                out_str += " "

            if newline:
                out_str += "\n"

            fid.write(out_str)

        # start with the banner
        self.__write_banner(fid, surf_name)
        fid.write(f"SURFACE\n")
        fid.write(f"{surf_name}\n")

        fid.write(f"#Nchordwise  Cspace  [Nspanwise  Sspace]\n")
        __write_data(["nchordwise", "cspace"], newline=False)
        if data['use surface spacing']:
            __write_data(["nspan", "sspace"])
        else:
            fid.write("\n")

        if "yduplicate" in data:
            fid.write("YDUPLICATE\n")
            __write_data(["yduplicate"])

        idx_surf = self.get_surface_index(surf_name)
        if idx_surf+1 != data["component"]:
            # only add component keys if we have to to avoid freaking 
            # people out who don't expect to see them
            fid.write("COMPONENT\n")
            __write_data(["component"])
        
        fid.write("SCALE\n")
        __write_data(["scale"])

        fid.write("TRANSLATE\n")
        __write_data(["translate"])

        fid.write("ANGLE\n")
        __write_data(["angle"])

        fid.write("#---------------------------------------\n")

        num_sec = data["chords"].size
        control_names = self.get_control_names()

        for idx_sec in range(num_sec):
            fid.write("SECTION\n")
            fid.write("#Xle      Yle      Zle      | Chord    Ainc     Nspan  Sspace\n")
            fid.write(
                f" {data['xyzles'][idx_sec, 0]:.6f} "
                f"{data['xyzles'][idx_sec, 1]:.6f} "
                f"{data['xyzles'][idx_sec, 2]:.6f}   "
                f"{data['chords'][idx_sec]:.6f} "
                f"{data['aincs'][idx_sec]:.6f} "
            )
            
            if data['nspans'][idx_sec] != 0:
                # if no section spacing is read in by avl nspans will be zero
                fid.write(
                    f"{data['nspans'][idx_sec]}      "
                    f"{data['sspaces'][idx_sec]}\n"
                )
            elif data['use surface spacing']:
                fid.write("\n")
            else:
                raise RuntimeError(f"neither surface nor section spacing information written"
                                   f"for surface {surf_name} section {idx_sec + 1}")

            afile = data["afiles"][idx_sec]

            if afile:
                fid.write(" AFILE\n")
                fid.write(f" {afile}\n")
            # output claf and  clcd if not default
            if data["claf"][idx_sec]!= 1.0:
                fid.write(" CLAF\n")
                fid.write(f" {data['claf'][idx_sec]}\n")

            if (data["clcdsec"][idx_sec] != 0.0).any():
                fid.write(" CDCL\n")
                fid.write(
                    f" {data['clcdsec'][idx_sec, 0]:.6f} "
                    f" {data['clcdsec'][idx_sec, 1]:.6f} "
                    f" {data['clcdsec'][idx_sec, 2]:.6f} "
                    f" {data['clcdsec'][idx_sec, 3]:.6f} "
                    f" {data['clcdsec'][idx_sec, 4]:.6f} "
                    f" {data['clcdsec'][idx_sec, 5]:.6f}\n"
                )
            

            # check for control surfaces

            for idx_local_cont_surf, idx_cont_surf in enumerate(data["icontd"][idx_sec]):
                fid.write(" CONTROL\n")
                fid.write("#surface   gain xhinge       hvec       SgnDup\n")
                fid.write(f" {control_names[idx_cont_surf-1]} ")
                fid.write(f" {data['gaind'][idx_sec][idx_local_cont_surf]}")
                fid.write(f" {data['xhinged'][idx_sec][idx_local_cont_surf]}")
                vhinge = data["vhinged"][idx_sec][idx_local_cont_surf]
                fid.write(f" {vhinge[0]:.6f} {vhinge[1]:.6f} {vhinge[2]:.6f}")
                fid.write(f" {data['refld'][idx_sec][idx_local_cont_surf]}\n")

    def __decodeFortranString(self, fort_string) -> str:
        # TODO: need a more general solution for |S<variable> type

        if fort_string.dtype == np.dtype("|S0"):
            # there are no characters in the sting to add
            return ""
        if fort_string.dtype == np.dtype("|S1"):
            py_string = b"".join(fort_string).decode().strip()
        elif fort_string.dtype == np.dtype("<U1"):
            py_string = "".join(fort_string).strip()
        elif fort_string.dtype == np.dtype("|S16"):
            py_string = fort_string.decode().strip()
        elif fort_string.dtype == np.dtype("|S40"):
            py_string = fort_string.decode().strip()
        elif fort_string.dtype == np.dtype("|S80"):
            py_string = fort_string.decode().strip()
        elif fort_string.dtype == np.dtype("|S120"):
            py_string = fort_string.decode().strip()
        else:
            raise TypeError(f"Unable to convert {fort_string} of type {fort_string.dtype} to string")

        return py_string

# region --- Utility functions
    def get_num_surfaces(self) -> int:
        """Get the number of surfaces in the geometry"""
        return self.get_avl_fort_arr("CASE_I", "NSURF")

    def get_surface_index(self, surf_name) -> int:
        surf_names = self.get_surface_names()
        idx_surf = surf_names.index(surf_name)
        return idx_surf

    def get_num_sections(self, surf_name) -> int:
        """Get the number of surfaces in the geometry"""
        idx_surf = self.get_surface_index(surf_name)
        slice_idx_surf = (idx_surf,)
        return self.get_avl_fort_arr("SURF_GEOM_I", "NSEC", slicer=slice_idx_surf)

    def get_num_strips(self) -> int:
        """Get the number of strips in the mesh"""
        return self.get_avl_fort_arr("CASE_I", "NSTRIP")

    def get_num_control_surfs(self) -> int:
        return self.get_avl_fort_arr("CASE_I", "NCONTROL")

    def get_mesh_size(self) -> int:
        """Get the number of panels in the mesh"""
        return int(self.get_avl_fort_arr("CASE_I", "NVOR"))

    def _createFortranStringArray(self, strList, num_max_char):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takesa list of strings and returns the array"""

        arr = np.zeros((len(strList), num_max_char), dtype="str")
        arr[:] = " "
        for i, s in enumerate(strList):
            for j in range(len(s)):
                arr[i, j] = s[j]

        return arr

    def _convertFortranStringArrayToList(self, fortArray):
        """Undoes the _createFotranStringArray"""
        strList = []

        if fortArray.size == 1:
            # we must handle the 0-d array case sperately
            return self.__decodeFortranString(fortArray[()])

        for ii in range(fortArray.size):
            py_string = self.__decodeFortranString(fortArray[ii])
            if py_string != "":
                strList.append(py_string)

        return strList
    
    def split_deriv_key(self, key):
        try:
            var, func = key.split('/')
        except Exception as e :
            print(key)
            import pdb; pdb.set_trace()
            raise e
        # remove leading 'd's
        var = var[1:]
        func = func[1:]
        return var, func
    
    def get_deriv_key(self, var: str, func : str) -> str:
        return f'd{func}/d{var}'
    
# ---------------------------
# --- Derivative routines ---
# ---------------------------

# --- input ad seeds --- 
    def get_constraint_ad_seeds(self) -> Dict[str, float]:
        con_seeds = {}
        for con in self.con_var_to_fort_var:
            idx_con = self.conval_idx_dict[con]
            blk = "CASE_R" + self.ad_suffix
            var = "CONVAL" + self.ad_suffix
            slicer = (0, idx_con)

            fort_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)
            con_seeds[con] = copy.deepcopy(fort_arr)

        return con_seeds

    def set_constraint_ad_seeds(self, con_seeds: Dict[str, Dict[str, float]], mode: str = "AD", scale=1.0) -> None:
        for con in con_seeds:
            # determine the proper index

            idx_con = self.conval_idx_dict[con]
            # determine the proper index

            con_seed_arr = con_seeds[con]

            if mode == "AD":
                # [0] because optvl only supports 1 run case
                blk = "CASE_R" + self.ad_suffix
                var = "CONVAL" + self.ad_suffix
                val = con_seed_arr * scale
                slicer = (0, idx_con)
                
                self.set_avl_fort_arr(blk, var, val, slicer=slicer)

            elif mode == "FD":
                # reverse lookup in the con_var_to_fort_var dict

                val = self.get_constraint(con)

                val += con_seed_arr * scale

                # use the contraint API to adjust the value
                self.set_constraint(con, val)
    
    def set_parameter_ad_seeds(self, parm_seeds: Dict[str, float], mode: str = "AD", scale=1.0) -> None:
        for param_key in parm_seeds:
            # blk, var = self.param_idx_dict[param_key]
            blk = "CASE_R"
            var = "PARVAL"
            
            idx_param = self.param_idx_dict[param_key]
            slicer = (0, idx_param)

            if mode == "AD":
                blk += self.ad_suffix
                var += self.ad_suffix
                val = parm_seeds[param_key] * scale
            elif mode == "FD":
                val = self.get_avl_fort_arr(blk, var, slicer=slicer)
                val += parm_seeds[param_key] * scale
            
            self.set_avl_fort_arr(blk, var, val, slicer=slicer)
 
    def get_parameter_ad_seeds(self) -> Dict[str, float]:
        # for param_key in parm_seeds:
            # blk, var = self.param_idx_dict[param_key]
        blk = "CASE_R"
        var = "PARVAL"
        
        slicer = (0, slice(None))

        blk += self.ad_suffix
        var += self.ad_suffix
        
        parval_seeds_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)
        
        param_seeds = {}
        for key, idx in self.param_idx_dict.items():
            param_seeds[key] = parval_seeds_arr[idx]
        
        return param_seeds
    
    def set_reference_ad_seeds(self, ref_seeds: Dict[str, float], mode: str = "AD", scale=1.0) -> None:
        for ref_key in ref_seeds:
            blk, var = self.ref_var_to_fort_var[ref_key]

            if mode == "AD":
                blk += self.ad_suffix
                var += self.ad_suffix
                val = ref_seeds[ref_key] * scale
            elif mode == "FD":
                val = self.get_avl_fort_arr(blk, var)
                val += ref_seeds[ref_key] * scale
            
            self.set_avl_fort_arr(blk, var, val)
            
    def get_reference_ad_seeds(self) -> Dict[str, float]:
        ref_seeds = {}
        for ref_key in self.ref_var_to_fort_var:
            blk, var = self.ref_var_to_fort_var[ref_key]
        
            blk += self.ad_suffix
            var += self.ad_suffix
            
            val = self.get_avl_fort_arr(blk, var)
            ref_seeds[ref_key] = copy.deepcopy(val)
        
        return ref_seeds
        
    def get_geom_ad_seeds(self) -> Dict[str, Dict[str, float]]:
        geom_seeds = {}
        for surf_key in self.surf_geom_to_fort_var:
            geom_seeds[surf_key] = {}
            for geom_key in self.surf_geom_to_fort_var[surf_key]:
                blk, var, slicer = self.surf_geom_to_fort_var[surf_key][geom_key]

                blk += self.ad_suffix
                var += self.ad_suffix

                geom_seeds[surf_key][geom_key] = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))

        return geom_seeds

    def set_geom_ad_seeds(self, geom_seeds: Dict[str, float], mode: str = "AD", scale=1.0) -> None:
        for surf_key in geom_seeds:
            for geom_key in geom_seeds[surf_key]:
                blk, var, slicer = self.surf_geom_to_fort_var[surf_key][geom_key]

                if mode == "AD":
                    blk += self.ad_suffix
                    var += self.ad_suffix
                    val = geom_seeds[surf_key][geom_key] * scale
                elif mode == "FD":
                    val = self.get_avl_fort_arr(blk, var, slicer=slicer)
                    val += geom_seeds[surf_key][geom_key] * scale
                # print(blk, var, val, slicer)
                self.set_avl_fort_arr(blk, var, val, slicer=slicer)

# --- state ad seeds ---
    def get_gamma_ad_seeds(self) -> np.ndarray:
        slicer = (slice(0, self.get_mesh_size()),)
        blk = "VRTX_R_DIFF"
        var = "GAM_DIFF"

        gamma_seeds = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))
        return gamma_seeds

    def set_gamma_ad_seeds(self, gamma_seeds: np.ndarray, mode: str = "AD", scale=1.0) -> None:
        slicer = (slice(0, self.get_mesh_size()),)
        if mode == "AD":
            blk = "VRTX_R_DIFF"
            var = "GAM_DIFF"
            val = gamma_seeds * scale
        elif mode == "FD":
            blk = "VRTX_R"
            var = "GAM"
            val = self.get_avl_fort_arr(blk, var, slicer=slicer)
            val += gamma_seeds * scale

        self.set_avl_fort_arr(blk, var, val, slicer=slicer)

    def get_gamma_d_ad_seeds(self) -> np.ndarray:
        slicer = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))
        blk = "VRTX_R_DIFF"
        var = "GAM_d_DIFF"

        gamma_d_seeds = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))
        return gamma_d_seeds

    def set_gamma_d_ad_seeds(self, gamma_d_seeds: np.ndarray, mode: str = "AD", scale=1.0) -> None:
        slicer = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))
        if mode == "AD":
            blk = "VRTX_R_DIFF"
            var = "GAM_D_DIFF"
            val = gamma_d_seeds * scale
        elif mode == "FD":
            blk = "VRTX_R"
            var = "GAM_D"
            val = self.get_avl_fort_arr(blk, var, slicer=slicer)
            val += gamma_d_seeds * scale

        self.set_avl_fort_arr(blk, var, val, slicer=slicer)

    def get_gamma_u_ad_seeds(self) -> np.ndarray:
        slicer = (slice(0, self.NUMAX), slice(0, self.get_mesh_size()))
        blk = "VRTX_R_DIFF"
        var = "GAM_U_DIFF"

        gamma_d_seeds = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))
        return gamma_d_seeds

    def set_gamma_u_ad_seeds(self, gamma_u_seeds: np.ndarray, mode: str = "AD", scale=1.0) -> None:
        slicer = (slice(0, self.NUMAX), slice(0, self.get_mesh_size()))
        if mode == "AD":
            blk = "VRTX_R_DIFF"
            var = "GAM_U_DIFF"
            val = gamma_u_seeds * scale
        elif mode == "FD":
            blk = "VRTX_R"
            var = "GAM_U"
            val = self.get_avl_fort_arr(blk, var, slicer=slicer)
            val += gamma_u_seeds * scale

        self.set_avl_fort_arr(blk, var, val, slicer=slicer)

# --- residual AD seeds --- 
    def get_residual_ad_seeds(self) -> np.ndarray:
        res_slice = (slice(0, self.get_mesh_size()),)
        res_seeds = copy.deepcopy(self.get_avl_fort_arr("VRTX_R_DIFF", "RES_DIFF", slicer=res_slice))
        return res_seeds

    def set_residual_ad_seeds(self, res_seeds: np.ndarray, scale=1.0) -> None:
        res_slice = (slice(0, self.get_mesh_size()),)
        self.set_avl_fort_arr("VRTX_R_DIFF", "RES_DIFF", res_seeds * scale, slicer=res_slice)
        return

    def get_residual_d_ad_seeds(self) -> np.ndarray:
        res_slice = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))
        res_seeds = copy.deepcopy(self.get_avl_fort_arr("VRTX_R_DIFF", "RES_D_DIFF", slicer=res_slice))
        return res_seeds

    def set_residual_d_ad_seeds(self, res_d_seeds: np.ndarray, scale=1.0) -> None:
        res_slice = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))

        self.set_avl_fort_arr("VRTX_R_DIFF", "RES_D_DIFF", res_d_seeds * scale, slicer=res_slice)
        return

    def get_residual_u_ad_seeds(self) -> np.ndarray:
        res_slice = (slice(0, self.NUMAX), slice(0, self.get_mesh_size()))
        res_seeds = copy.deepcopy(self.get_avl_fort_arr("VRTX_R_DIFF", "RES_U_DIFF", slicer=res_slice))
        return res_seeds

    def set_residual_u_ad_seeds(self, res_u_seeds: np.ndarray, scale=1.0) -> None:
        res_u_slice = (slice(0, self.NUMAX), slice(0, self.get_mesh_size()))

        self.set_avl_fort_arr("VRTX_R_DIFF", "RES_U_DIFF", res_u_seeds * scale, slicer=res_u_slice)
        return

# --- output AD seeds ---
    def get_function_ad_seeds(self):
        func_seeds = {}
        for _var in self.case_var_to_fort_var:
            blk, var = self.case_var_to_fort_var[_var]
            blk += self.ad_suffix
            var += self.ad_suffix
            val = self.get_avl_fort_arr(blk, var)
            func_seeds[_var] = copy.deepcopy(val)

        return func_seeds

    def set_function_ad_seeds(self, func_seeds: Dict[str, float], scale=1.0):
        for _var in func_seeds:
            blk, var = self.case_var_to_fort_var[_var]
            blk += self.ad_suffix
            var += self.ad_suffix
            val = func_seeds[_var] * scale
            self.set_avl_fort_arr(blk, var, val)

    def get_consurf_derivs_ad_seeds(self):
        cs_deriv_seeds = {}
        consurf_names = self.get_control_names()
        num_consurf = len(consurf_names)

        for _var in self.case_derivs_to_fort_var:
            blk, var = self.case_derivs_to_fort_var[_var]
            slicer = (slice(0, num_consurf),)
            blk += self.ad_suffix
            var += self.ad_suffix
            val_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)

            for idx_control, val in enumerate(val_arr):
                control = consurf_names[idx_control]
                cs_deriv_seeds[self.get_deriv_key(control,_var)] = val[()]

        return cs_deriv_seeds

    def set_consurf_derivs_ad_seeds(self, cs_deriv_seeds: Dict[str, float], scale=1.0):
        consurf_names = self.get_control_names()
        num_consurf = len(consurf_names)

        for deriv_func in cs_deriv_seeds:
            val_arr = np.zeros((num_consurf))
            
            var, cs_name = self.split_deriv_key(deriv_func)
            
            idx_cs = consurf_names.index(cs_name)
            val_arr[idx_cs] = cs_deriv_seeds[deriv_func] * scale

            blk, var = self.case_derivs_to_fort_var[var]
            slicer = (slice(0, num_consurf),)

            blk += self.ad_suffix
            var += self.ad_suffix

            self.set_avl_fort_arr(blk, var, val_arr, slicer=slicer)

    def get_stab_derivs_ad_seeds(self):
        
        deriv_data = {}

        for func_key, avl_vars in self.case_stab_derivs_to_fort_var.items():
            deriv_data[func_key] = {}

            blk, var = avl_vars
            blk += self.ad_suffix
            var += self.ad_suffix    
            val_arr = self.get_avl_fort_arr(blk, var)
            deriv_data[func_key] = val_arr[()]

        return deriv_data

    def set_stab_derivs_ad_seeds(self, stab_deriv_seeds: Dict[str, Dict[str, float]], scale=1.0):
        for func_key in stab_deriv_seeds:
            blk, var  = self.case_stab_derivs_to_fort_var[func_key]

            blk += self.ad_suffix
            var += self.ad_suffix
            
            val = stab_deriv_seeds[func_key] * scale
            
            self.set_avl_fort_arr(blk, var, val)

# --- derivative utils
    def clear_ad_seeds(self):
        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    if not (_var.startswith("__") and _var.endswith("__")):
                        val = getattr(diff_blk, _var)
                        setattr(diff_blk, _var, val * 0.0)

    def clear_ad_seeds_fast(self):
        # Only clear the seeds that are used in Make_tapenade file
        num_vor = self.get_mesh_size()
        gam = self.get_avl_fort_arr("VRTX_R", "GAM")
        num_vor_max = gam.size

        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    if not (_var.startswith("__") and _var.endswith("__")):

                        val = getattr(diff_blk, _var)

                        # trim sizes set to NVMAX to NVOR
                        shape = val.shape
                        slices = []
                        for idx_dim in range(len(shape)):
                            dim_size = shape[idx_dim]
                            if dim_size == num_vor_max:
                                dim_size = num_vor

                            slices.append(slice(0, dim_size))
                        slicer = tuple(slices)
                        val[slicer] = 0.0
                        # setattr(diff_blk, _var, val)
                        # print(diff_blk, _var, val.shape, slicer)
                        # mb_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
                        # print(f"    Memory usage: {mb_memory} MB")

    def print_ad_seeds(self, print_non_zero: bool = False):
        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    if not (_var.startswith("__") and _var.endswith("__")):
                        val = getattr(diff_blk, _var)
                        norm = np.linalg.norm(val)
                        if norm == 0.0 and print_non_zero:
                            continue

                        print(att, _var, norm)

# --- jacobian vecotr products --- 
    def execute_jac_vec_prod_fwd(
        self,
        con_seeds: Optional[Dict[str, float]] = None,
        geom_seeds: Optional[Dict[str, Dict[str, any]]] = None,
        param_seeds: Optional[Dict[str, float]] = None,
        ref_seeds: Optional[Dict[str, float]] = None,
        gamma_seeds: Optional[np.ndarray] = None,
        gamma_d_seeds: Optional[np.ndarray] = None,
        gamma_u_seeds: Optional[np.ndarray] = None,
        mode="AD",
        step=1e-7,
    ) -> Tuple[Dict[str, float], np.ndarray]:
        # TODO: add better name for gamma_d, it is too confusing

        # TODO: add error if data is not initailzed properly
        #   The easiest fix is to run an analysis before hand

        mesh_size = self.get_mesh_size()
        num_control_surfs = self.get_num_control_surfs()

        if con_seeds is None:
            con_seeds = {}

        if geom_seeds is None:
            geom_seeds = {}

        if gamma_seeds is None:
            gamma_seeds = np.zeros(mesh_size)

        if gamma_d_seeds is None:
            gamma_d_seeds = np.zeros((num_control_surfs, mesh_size))
        
        if gamma_u_seeds is None:
            gamma_u_seeds = np.zeros((self.NUMAX, mesh_size))
        
        if param_seeds is None:
            param_seeds = {}
            
        if ref_seeds is None:
            ref_seeds = {}
        
        res_slice = (slice(0, mesh_size),)
        res_d_slice = (slice(0, num_control_surfs), slice(0, mesh_size))
        res_u_slice = (slice(0, self.NUMAX), slice(0, mesh_size))

        if mode == "AD":
            # set derivative seeds
            # self.clear_ad_seeds()
            self.set_constraint_ad_seeds(con_seeds)
            self.set_geom_ad_seeds(geom_seeds)
            self.set_gamma_ad_seeds(gamma_seeds)
            self.set_gamma_d_ad_seeds(gamma_d_seeds)
            self.set_gamma_u_ad_seeds(gamma_u_seeds)
            self.set_parameter_ad_seeds(param_seeds)
            self.set_reference_ad_seeds(ref_seeds)

            self.avl.update_surfaces_d()
            self.avl.get_res_d()
            self.avl.velsum_d()
            self.avl.aero_d()

            # extract derivatives seeds and set the output dict of functions
            func_seeds = self.get_function_ad_seeds()
            res_seeds = self.get_residual_ad_seeds()
            consurf_derivs_seeds = self.get_consurf_derivs_ad_seeds()
            stab_derivs_seeds = self.get_stab_derivs_ad_seeds()
            res_d_seeds = self.get_residual_d_ad_seeds()
            res_u_seeds = self.get_residual_u_ad_seeds()

            self.set_constraint_ad_seeds(con_seeds, scale=0.0)
            self.set_geom_ad_seeds(geom_seeds, scale=0.0)
            self.set_gamma_ad_seeds(gamma_seeds, scale=0.0)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, scale=0.0)
            self.set_gamma_u_ad_seeds(gamma_u_seeds, scale=0.0)
            self.set_parameter_ad_seeds(param_seeds, scale=0.0)
            self.set_reference_ad_seeds(ref_seeds, scale=0.0)

            #TODO: remove??
            self.set_avl_fort_arr("VRTX_R_DIFF", "GAM_DIFF", gamma_seeds * 0.0, slicer=res_slice)

        if mode == "FD":
            self.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
            self.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)
            self.set_gamma_ad_seeds(gamma_seeds, mode="FD", scale=step)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, mode="FD", scale=step)
            self.set_gamma_u_ad_seeds(gamma_u_seeds, mode="FD", scale=step)
            self.set_parameter_ad_seeds(param_seeds, mode="FD", scale=step)
            self.set_reference_ad_seeds(ref_seeds, mode="FD", scale=step)

            # propogate the seeds through without resolving
            self.avl.update_surfaces()
            self.avl.get_res()
            self.avl.velsum()
            self.avl.aero()

            coef_data_peturb = self.get_total_forces()
            consurf_derivs_petrub = self.get_control_stab_derivs()
            stab_deriv_petrub = self.get_stab_derivs()

            res_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES", slicer=res_slice))
            res_d_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_D", slicer=res_d_slice))
            res_u_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_U", slicer=res_u_slice))
            
            
            self.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-1 * step)
            self.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_ad_seeds(gamma_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_u_ad_seeds(gamma_u_seeds, mode="FD", scale=-1 * step)
            self.set_parameter_ad_seeds(param_seeds, mode="FD", scale=-1*step)
            self.set_reference_ad_seeds(ref_seeds, mode="FD", scale=-1*step)

            self.avl.update_surfaces()
            self.avl.get_res()
            self.avl.velsum()
            self.avl.aero()

            coef_data = self.get_total_forces()
            consurf_derivs = self.get_control_stab_derivs()
            stab_deriv = self.get_stab_derivs()

            res = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES", slicer=res_slice))
            res_d = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_D", slicer=res_d_slice))
            res_u = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_U", slicer=res_u_slice))

            func_seeds = {}
            for func_key in coef_data:
                func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

            consurf_derivs_seeds = {}
            for deriv_func in consurf_derivs:
                consurf_derivs_seeds[deriv_func] = (
                    consurf_derivs_petrub[deriv_func] - consurf_derivs[deriv_func]
                ) / step
                    
            stab_derivs_seeds = {}
            for deriv_func in stab_deriv:
                stab_derivs_seeds[deriv_func] = (
                    stab_deriv_petrub[deriv_func] - stab_deriv[deriv_func]
                ) / step

            res_seeds = (res_peturbed - res) / step
            res_d_seeds = (res_d_peturbed - res_d) / step
            res_u_seeds = (res_u_peturbed - res_u) / step
        
        # TODO-clean: the way these arrays are returned is a bit of a mess 
        return func_seeds, res_seeds, consurf_derivs_seeds, stab_derivs_seeds, res_d_seeds, res_u_seeds

    def execute_jac_vec_prod_rev(
        self,
        func_seeds: Optional[Dict[str, float]] = None,
        res_seeds: Optional[np.ndarray] = None,
        consurf_derivs_seeds: Optional[Dict[str, float]] = None,
        stab_derivs_seeds: Optional[Dict[str, float]] = None,
        res_d_seeds: Optional[np.ndarray] = None,
        res_u_seeds: Optional[np.ndarray] = None,
        print_timings=False,
    ) -> Tuple[Dict[str, float], Dict[str, Dict[str, any]], np.ndarray]:
        # extract derivatives seeds and set the output dict of functions
        mesh_size = self.get_mesh_size()
        num_surf = self.get_num_control_surfs()

        if func_seeds is None:
            func_seeds = {}

        if res_seeds is None:
            res_seeds = np.zeros(mesh_size)

        if res_d_seeds is None:
            res_d_seeds = np.zeros((num_surf, mesh_size))
        
        if res_u_seeds is None:
            res_u_seeds = np.zeros((self.NUMAX, mesh_size))

        if consurf_derivs_seeds is None:
            consurf_derivs_seeds = {}
        
        if stab_derivs_seeds is None:
            stab_derivs_seeds = {}

        # set derivative seeds
        # self.clear_ad_seeds()
        time_last = time.time()
        self.set_function_ad_seeds(func_seeds)
        self.set_residual_ad_seeds(res_seeds)
        self.set_residual_d_ad_seeds(res_d_seeds)
        self.set_residual_u_ad_seeds(res_u_seeds)
        self.set_consurf_derivs_ad_seeds(consurf_derivs_seeds)
        self.set_stab_derivs_ad_seeds(stab_derivs_seeds)

        if print_timings:
            print(f"    Time to set seeds: {time.time() - time_last}")
            time_last = time.time()

        # propogate the seeds through without resolveing
        self.avl.aero_b()
        self.avl.velsum_b()
        self.avl.get_res_b()
        self.avl.update_surfaces_b()
        if print_timings:
            print(f"    Time to propogate seeds: {time.time() - time_last}")
            time_last = time.time()

        # extract derivatives seeds and set the output dict of functions
        con_seeds = self.get_constraint_ad_seeds()
        geom_seeds = self.get_geom_ad_seeds()
        gamma_seeds = self.get_gamma_ad_seeds()
        gamma_d_seeds = self.get_gamma_d_ad_seeds()
        gamma_u_seeds = self.get_gamma_u_ad_seeds()
        param_seeds = self.get_parameter_ad_seeds()
        ref_seeds = self.get_reference_ad_seeds()
        if print_timings:
            print(f"    Time to extract seeds: {time.time() - time_last}")
            time_last = time.time()

        self.set_function_ad_seeds(func_seeds, scale=0.0)
        self.set_residual_ad_seeds(res_seeds, scale=0.0)
        self.set_residual_d_ad_seeds(res_d_seeds, scale=0.0)
        self.set_residual_u_ad_seeds(res_u_seeds, scale=0.0)
        self.set_consurf_derivs_ad_seeds(consurf_derivs_seeds, scale=0.0)
        self.set_stab_derivs_ad_seeds(stab_derivs_seeds, scale=0.0)
        if print_timings:
            print(f"    Time to clear seeds: {time.time() - time_last}")
            time_last = time.time()

        return con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds, gamma_u_seeds, param_seeds, ref_seeds

    def execute_adjoint_solve():
        raise NotImplementedError

    def execute_direct_solve():
        raise NotImplementedError

    def execute_run_sensitivies(self, funcs, stab_derivs=None, consurf_derivs=None, print_timings=False):

        """
        only runs in adjoint mode
        funcs: list of functions that need derivatives
        """
        sens = {}
        
        if self.get_avl_fort_arr("CASE_L", "LTIMING"):
            print_timings = True
        
        # set up and solve the adjoint for each function
        for func in funcs:
            sens[func] = {}
            # get the RHS of the adjoint equation (pFpU)
            # TODO: remove seeds if it doesn't effect accuracy
            # self.clear_ad_seeds()
            time_last = time.time()
            _, _, pfpU, _, _, _, _ = self.execute_jac_vec_prod_rev(func_seeds={func: 1.0})
            if print_timings:
                print(f"Time to get RHS: {time.time() - time_last}")
                time_last = time.time()

            # self.clear_ad_seeds()
            # u solver adjoint equation with RHS
            self.set_gamma_ad_seeds(-1 * pfpU)
            solve_stab_deriv_adj=False
            solve_con_surf_adj=False
            self.avl.solve_adjoint(solve_stab_deriv_adj, solve_con_surf_adj)
            if print_timings:
                print(f"Time to solve adjoint: {time.time() - time_last}")
                time_last = time.time()
            # get the resulting adjoint vector (dfunc/dRes) from fortran
            dfdR = self.get_residual_ad_seeds()
            # self.clear_ad_seeds()
            con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self.execute_jac_vec_prod_rev(
                func_seeds={func: 1.0}, res_seeds=dfdR
            )
            if print_timings:
                print(f"Time to combine derivs: {time.time() - time_last}")
                time_last = time.time()

            sens[func].update(con_seeds)
            sens[func].update(geom_seeds)
            sens[func].update(param_seeds)
            sens[func].update(ref_seeds)

        if consurf_derivs is not None:
            if print_timings:
                print("Running consurf derivs")
                time_last = time.time()

            for func_key in consurf_derivs:
                if func_key not in sens:
                    sens[func_key] = {}

                # get the RHS of the adjoint equation (pFpU)
                # TODO: remove seeds if it doesn't effect accuracy
                _, _, pfpU, pf_pU_d, _, _, _ = self.execute_jac_vec_prod_rev(consurf_derivs_seeds={func_key: 1.0})
                if print_timings:
                    print(f"Time to get RHS: {time.time() - time_last}")
                    time_last = time.time()

                # self.clear_ad_seeds()
                # u solver adjoint equation with RHS
                self.set_gamma_ad_seeds(-1 * pfpU)
                self.set_gamma_d_ad_seeds(-1 * pf_pU_d)
                solve_stab_deriv_adj=False
                solve_con_surf_adj=True
                self.avl.solve_adjoint(solve_stab_deriv_adj, solve_con_surf_adj)
                if print_timings:
                    print(f"Time to solve adjoint: {time.time() - time_last}")
                    time_last = time.time()

                # get the resulting adjoint vector (dfunc/dRes) from fortran
                dfdR = self.get_residual_ad_seeds()
                dfdR_d = self.get_residual_d_ad_seeds()
                # self.clear_ad_seeds()
                con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self.execute_jac_vec_prod_rev(
                    consurf_derivs_seeds={func_key: 1.0}, res_seeds=dfdR, res_d_seeds=dfdR_d
                )
                if print_timings:
                    print(f"Time to combine : {time.time() - time_last}")
                    time_last = time.time()

                sens[func_key].update(con_seeds)
                sens[func_key].update(geom_seeds)
                sens[func_key].update(param_seeds)
                sens[func_key].update(ref_seeds)

        if stab_derivs is not None:
            if print_timings:
                print("Running stab derivs")
                time_last = time.time()

            # sd_deriv_seeds = {}
            for func_key in stab_derivs:
                # sd_deriv_seeds[func_key] = {}
                if func_key not in sens:
                    sens[func_key] = {}

                # for var_key in stab_derivs[func_key]:
                #     sd_deriv_seeds[func_key][var_key] = 1.0
                #     sens[func_key][var_key] = {}

                # get the RHS of the adjoint equation (pFpU)
                # TODO: remove seeds if it doesn't effect accuracy
                _, _, pfpU, _, pf_pU_u, _, _ = self.execute_jac_vec_prod_rev(stab_derivs_seeds={func_key : 1.0})
                if print_timings:
                    print(f"Time to get RHS: {time.time() - time_last}")
                    time_last = time.time()

                # self.clear_ad_seeds()
                # u solver adjoint equation with RHS
                self.set_gamma_ad_seeds(-1 * pfpU)
                self.set_gamma_u_ad_seeds(-1 * pf_pU_u)
                solve_stab_deriv_adj=True
                solve_con_surf_adj=False
                self.avl.solve_adjoint(solve_stab_deriv_adj, solve_con_surf_adj)
                if print_timings:
                    print(f"Time to solve adjoint: {time.time() - time_last}")
                    time_last = time.time()

                # get the resulting adjoint vector (dfunc/dRes) from fortran
                dfdR = self.get_residual_ad_seeds()
                dfdR_u = self.get_residual_u_ad_seeds()
                # self.clear_ad_seeds()
                con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self.execute_jac_vec_prod_rev(
                    stab_derivs_seeds={func_key : 1.0}, res_seeds=dfdR, res_u_seeds=dfdR_u
                )

                if print_timings:
                    print(f"Time to combine : {time.time() - time_last}")
                    time_last = time.time()

                sens[func_key].update(con_seeds)
                sens[func_key].update(geom_seeds)
                sens[func_key].update(param_seeds)
                sens[func_key].update(ref_seeds)
                # sd_deriv_seeds[func_key] = 0.0

        return sens

# --- ploting and vizulaization ---
    def add_mesh_plot(self, axis, xaxis='x', yaxis='y', show_mesh=True):
        """ adds a plot of the aircraft mesh to the axis"""
        mesh_size = self.get_mesh_size()
        num_control_surfs = self.get_num_control_surfs()
        num_strips = self.get_num_strips()
        num_surfs = self.get_num_surfaces()
            
        # get the mesh points for ploting
        mesh_slice = (slice(0, mesh_size),)
        strip_slice = (slice(0, num_strips),)
        surf_slice = (slice(0, num_surfs),)
        
        rv1 = self.get_avl_fort_arr("VRTX_R", "RV1", slicer=mesh_slice)
        rv2 = self.get_avl_fort_arr("VRTX_R", "RV2", slicer=mesh_slice)
        rle1 = self.get_avl_fort_arr("STRP_R", "RLE1", slicer=strip_slice)
        rle2 = self.get_avl_fort_arr("STRP_R", "RLE2", slicer=strip_slice)
        chord1 = self.get_avl_fort_arr("STRP_R", "CHORD1", slicer=strip_slice)
        chord2 = self.get_avl_fort_arr("STRP_R", "CHORD2", slicer=strip_slice)
        jfrst = self.get_avl_fort_arr("SURF_I", "JFRST", slicer=surf_slice)
        
        ijfrst = self.get_avl_fort_arr("STRP_I", "IJFRST", slicer=strip_slice)
        nvstrp = self.get_avl_fort_arr("STRP_I", "NVSTRP", slicer=strip_slice)
        
        
        nj = self.get_avl_fort_arr("SURF_I", "NJ", slicer=surf_slice)
        imags = self.get_avl_fort_arr("SURF_I", "IMAGS")
              
        for idx_surf in range(num_surfs):
            # get the range of the elements that belong to this surfaces
            strip_st = jfrst[idx_surf] -1
            strip_end = strip_st + nj[idx_surf]
            
            # inboard and outboard of outline
            # get surfaces that have not been duplicated
            if imags[idx_surf] > 0:
                j1 = strip_st
                jn = strip_end  - 1
                dj=1
            else: 
                # this surface is a duplicate
                j1 = strip_end  - 1
                jn = strip_st
                dj = -1
                
            pts = {
                'x':[rle1[j1, 0], rle1[j1, 0] + chord1[j1]],
                'y': [rle1[j1, 1], rle1[j1, 1]],
                'z': [rle1[j1, 2], rle1[j1, 2]]         
            }
            # # chord-wise grid
            axis.plot(pts[xaxis], pts[yaxis], color='black')
            
                
            pts = {
                'x' :np.array([rle2[jn, 0], rle2[jn, 0] + chord2[jn]]),
                'y' :np.array([rle2[jn, 1], rle2[jn, 1]]),
                'z' :np.array([rle2[jn, 2], rle2[jn, 2]]),
            }
            
            # # chord-wise grid
            axis.plot(pts[xaxis], pts[yaxis], color='black')
            
            # # --- outline of surface ---
            # front 
            pts ={
                'x': np.append(rle1[j1:jn:dj, 0], rle2[jn, 0]),
                'y': np.append(rle1[j1:jn:dj, 1], rle2[jn, 1]),
                'z': np.append(rle1[j1:jn:dj, 2], rle2[jn, 2])
            }
            axis.plot(pts[xaxis], pts[yaxis], '-', color='black')
            
            # aft
            
            pts = {
                'x':  np.append(rle1[j1:jn:dj, 0]+ chord1[j1:jn:dj], rle2[jn, 0]+chord2[jn]),
                'y':  np.append(rle1[j1:jn:dj, 1], rle2[jn, 1]),
                'z':  np.append(rle1[j1:jn:dj, 2], rle2[jn, 2]),
            }
            axis.plot(pts[xaxis], pts[yaxis], '-', color='black')
            
            
            if show_mesh:
                for idx_strip in range(strip_st, strip_end):
                    if idx_strip != strip_st:
                        
                        pts = {
                            'x':[rle1[idx_strip, 0], rle1[idx_strip, 0] + chord1[idx_strip]],
                            'y':[rle1[idx_strip, 1], rle1[idx_strip, 1]],
                            'z':[rle1[idx_strip, 2], rle1[idx_strip, 2]],
                        }         
                        
                        # # chord-wise grid
                        axis.plot(pts[xaxis], pts[yaxis], '--', color='grey', linewidth=0.3)
                        
                    
                    
                    vor_st = ijfrst[idx_strip] -1
                    vor_end = vor_st + nvstrp[idx_strip]

                    # spanwise grid
                    for idx_vor in range(vor_st,vor_end):
                        pts = {
                            'x': [rv1[idx_vor, 0], rv2[idx_vor, 0]],
                            'y': [rv1[idx_vor, 1], rv2[idx_vor, 1]],
                            'z': [rv1[idx_vor, 2], rv2[idx_vor, 2]],
                        }
                        axis.plot(pts[xaxis], pts[yaxis], '--', color='grey', linewidth=0.3)
                        
    def plot_geom(self, **kwargs):
        import matplotlib.pyplot as plt
        
        ax = plt.subplot(2,1,1)
        self.add_mesh_plot(ax, xaxis='y', yaxis='x')
        ax.set_ylabel('X', rotation=0)
        ax = plt.subplot(2,1,2)
        self.add_mesh_plot(ax, xaxis='y', yaxis='z')
        ax.set_ylabel('Z', rotation=0)
        ax.set_xlabel('Y')
    
        plt.axis('equal')
        plt.show()        

    def get_cp_data(self):
        nvort = self.get_mesh_size()
        mesh_slicer = (slice(0, nvort),)
        
        num_surfs = self.get_num_surfaces()
        surf_slice = (slice(0,num_surfs))
        
        nj = self.get_avl_fort_arr('SURF_I', 'NJ', slicer=surf_slice)
        nk = self.get_avl_fort_arr('SURF_I', 'NK', slicer=surf_slice)
        
        
        xyz_list = []
        cp_list = []
        
        for idx_surf in range(num_surfs):
            nChords = nk[idx_surf]
            nStrips = nj[idx_surf]
            nPts = (nStrips+1)*(nChords*2+1)
            nCCPts = (nStrips)*(nChords*2)
            
            vtx_slicer = (idx_surf, slice(0, nPts),slice(None))
            cp_slicer = (idx_surf, slice(0, nCCPts))


            xyz = self.get_avl_fort_arr('VRTX_S', 'XYZSURF', slicer=vtx_slicer)
            
            xyz = xyz.reshape((nStrips+1, nChords*2+1, 3))
            xyz_list.append(xyz)
        
            cp = self.get_avl_fort_arr('VRTX_S', 'CPSURF', slicer=cp_slicer)
            cp = cp.reshape((nStrips, nChords*2))
            cp_list.append(cp)
        
        return xyz_list, cp_list

    def plot_cp(self):
        import matplotlib.pyplot as plt
        from matplotlib import cm
        
        self.avl.cpoml(False)
        xyz, cp = self.get_cp_data()
        
        num_surfs = self.get_num_surfaces()
        
        # create the map for the cp color
        cp_max = -1e99
        cp_min = 1e99
        for idx_surf in range(num_surfs):
            cp_max = max(cp_max, np.max(cp[idx_surf]))
            cp_min = min(cp_min, np.min(cp[idx_surf]))

        cp_amax = max(np.abs(cp_min), np.abs(cp_max))
        # Create a normalized colormap
        norm = plt.Normalize(vmin=-1*cp_amax, vmax=cp_amax)
        m = cm.ScalarMappable(cmap=cm.bwr)
        m.set_clim(vmin=cp_min, vmax=cp_max)
        
        # do the actual ploting of each surface 
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        plt.subplots_adjust(left=0.025, right=0.925, top=0.925, bottom=0.025)

        
        for idx_surf in range(num_surfs):
                    
            xyz_surf = xyz[idx_surf]
            
            face_color = m.to_rgba(cp[idx_surf])  
            
            ax.plot_surface(xyz_surf[:,:,0],xyz_surf[:,:,1], xyz_surf[:,:,2], facecolors = face_color)
        
        
        plt.axis('off')

        plt.grid(b=None)
        fig.colorbar(m, ax=ax)
        # Set an equal aspect ratio
        ax.set_aspect('equal')

        plt.show()

    def write_tecplot(self, file_name : str, solution_time=None):
        if solution_time is not None:
            add_time = True
        else:
            solution_time = 0.0
            add_time = False
            
        self.avl.cpoml(False)
        self.avl.write_tecplot(file_name + '.dat', add_time, solution_time )
        
