"""
This module contains the OVLSolver, which is the main wrapper for the Fortran level API.
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import copy
from typing import Dict, List, Tuple, Any, TextIO
import warnings
import glob
from typing import Optional
import platform
from collections import OrderedDict
import operator
from itertools import chain
import importlib.metadata

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np


# =============================================================================
# Extension modules
# =============================================================================
from . import MExt
from .utils.check_surface_dict import pre_check_input_dict
from .utils.airfoil_utils import read_coordinates_files

# import the version of the package to include in some io and help messages
__version__ = importlib.metadata.version(__package__ or __name__)

class OVLSolver(object):
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
        "visc Cm_a": 28,
        "visc Cm_u": 29,
    }

    # fmt: off
    # This dict has the following structure:
    # python key: [common block name, fortran varaiable name]
    case_var_to_fort_var = {
        # lift and drag from surface integration (wind frame)
        "CL": ["CASE_R", "CLTOT"],
        "CD": ["CASE_R", "CDTOT"],
        "CDv": ["CASE_R", "CDVTOT"], # viscous drag
        "CDi": ["CASE_R", "CDITOT"], #  induced drag
        
        # lift and drag calculated from farfield integration
        "CLff": ["CASE_R", "CLFF"],
        "CYff": ["CASE_R", "CYFF"],
        "CDff": ["CASE_R", "CDFF"],
        
        # # non-dimensionalized forces 
        # "CX": ["CASE_R", "CFTOT", 0], 
        # "CY": ["CASE_R", "CFTOT", 1], 
        # "CZ": ["CASE_R", "CFTOT", 2],   
        
        # non-dimensionalized forces (body frame)
        "CX": ["CASE_R", "CXBAX"], 
        "CY": ["CASE_R", "CYBAX"], 
        "CZ": ["CASE_R", "CZBAX"],   
        
        # non-dimensionalized moments (body frame)
        "Cl": ["CASE_R", "CRBAX"],
        "Cm": ["CASE_R", "CMBAX"],
        "Cn": ["CASE_R", "CNBAX"],

        # non-dimensionalized moments (stablity frame)
        "Cl'": ["CASE_R", "CRSAX"],
        "Cn'": ["CASE_R", "CNSAX"],
        
        # spanwise efficiency
        "e": ["CASE_R", "SPANEF"],
    }
    
    body_forces_var_to_fort_var = {
        'length':        ['BODY_R', 'ELBDY'],
        'surface area':  ['BODY_R', 'SRFBDY'],
        'volume':        ['BODY_R', 'VOLBDY'],
        'CL':            ['BODY_R', 'CLBDY'],
        'CD':            ['BODY_R', 'CDBDY'],
        'Cm':            ['BODY_R', 'CMBDY', 1],
        'CY':            ['BODY_R', 'CFBDY', 1],
        'Cn':            ['BODY_R', 'CMBDYBAX', 2],
        'Cl':            ['BODY_R', 'CMBDYBAX', 0]
    }
    
    ref_var_to_fort_var = {
        "Sref": ["CASE_R", "SREF"],
        "Cref": ["CASE_R", "CREF"],
        "Bref": ["CASE_R", "BREF"],
        "XYZref": ["CASE_R", "XYZREF"],
    }

    case_derivs_to_fort_var = {
        # derivative of coefficents wrt control surface deflections
        
        # stablity axis values
        "CL": ["CASE_R", "CLTOT_D"],
        "CD": ["CASE_R", "CDTOT_D"],
        "Cl'": ["CASE_R", "CRSAX_D"],
        # "Cm'" is the same as "Cm"
        "Cn'": ["CASE_R", "CNSAX_D"],
        
        # body axis values
        "CX": ["CASE_R", "CXTOT_D_BA"],
        "CY": ["CASE_R", "CYTOT_D_BA"],
        "CZ": ["CASE_R", "CZTOT_D_BA"],
        "Cl": ["CASE_R", "CRTOT_D_BA"],
        "Cm": ["CASE_R", "CMTOT_D_BA"],
        "Cn": ["CASE_R", "CNTOT_D_BA"],
    }

    # This dict has the following structure:
    # python key: [common block name, fortran varaiable name]
    case_surf_var_to_fort_var = {
        # surface contributions to total lift and drag from surface integration (wind frame)
        "CL": ["SURF_R", "CLSURF"],
        "CD": ["SURF_R", "CDSURF"],
        "CDv": ["SURF_R", "CDVSURF"],  # viscous drag
        "CDi": ["SURF_R", "CDISURF"],  # inviscid drag
        
        # geometric
        "area": ["SURF_R", "SSURF"], 
        "average chord": ["SURF_R", "CAVESURF"],
        
        # non-dimensionalized forces 
        "CX": ["SURF_R", "CXSURF"], 
        "CY": ["SURF_R", "CYSURF"], 
        "CZ": ["SURF_R", "CZSURF"],   
         
        # non-dimensionalized momen ts (body frame)
        "Cl": ["SURF_R", "CRSURF"  ],
        "Cm": ["SURF_R", "CMSUR  F"],
        "Cn": ["SURF_R", "CNSURF"],
        
        # forces non-dimentionalized by surface quantities
        # uses surface area instead of sref and takes moments about leading edge
        "CL surf" : ["SURF_R", "CL_LSRF"],
        "CD surf" : ["SURF_R", "CD_LSRF"],
        "CMLE_LSTRP surf" : ["SURF_R", "CMLE_LSRF"],
        
        #TODO: add CF_LSRF(3,NFMAX), CM_LSRF(3,NFMAX)
    }

    body_geom_to_fort_var = {
        "scale": ["BODY_GEOM_R", "XYZSCAL_B"],
        "translate": ["BODY_GEOM_R", "XYZTRAN_B"],
        "yduplicate": ["BODY_GEOM_R", "YDUPL_B"],
        "bfile": ["CASE_C", "BFILES"],
        "nvb": ["BODY_GEOM_I", "NVB"],
        "bspace": ["BODY_GEOM_R", "BSPACE"],
    }


    header_var_to_fort_var = {
        "title": ["CASE_C", "TITLE"],
        "mach":  ["CASE_R", "MACH0"],
        "iysym": ["CASE_I", "IYSYM"],
        "izsym": ["CASE_I", "IZSYM"],
        "zsym":  ["CASE_R", "ZSYM"],
        "Sref": ["CASE_R", "SREF"],
        "Cref": ["CASE_R", "CREF"],
        "Bref": ["CASE_R", "BREF"],
        "XYZref":  ["CASE_R","XYZREF0"],
        "CDp":   ["CASE_R", "CDREF0"],
    }

    # fmt: on

    ad_suffix = "_DIFF"

    # Primary array limits: These also need to updated in the Fortran layer if changed
    NSMAX = 400  # number of chord strips
    NFMAX = 30  # number of surfaces
    NLMAX = 500  # number of source/doublet line nodes
    NBMAX = 20  # number of bodies
    NUMAX = 6  # number of freestream parameters (V,Omega)
    NDMAX = 30  # number of control deflection parameters
    NGMAX = 20  # number of design variables
    NRMAX = 25  # number of stored run cases
    NTMAX = 5000  # number of stored time levels
    IBX = 300
    ICONX = 20

    if platform.system == "Windows":
        NVMAX = 5000  # number of horseshoe vortices
    else:
        NVMAX = 6000  # number of horseshoe vortices

    def __init__(
        self,
        geo_file: Optional[str] = None,
        mass_file: Optional[str] = None,
        input_dict: Optional[dict] = None,
        debug: Optional[bool] = False,
        timing: Optional[bool] = False,
    ):
        """Initalize the python and fortran libary from the given objects

        Args:
            geo_file: AVL geometry file
            mass_file: AVL mass file
            debug: flag for debug printing
            timing: flag for timing printing

        """

        if timing:
            start_time = time.time()

        self.debug = debug

        # MExt is important for creating multiple instances of the AVL solver that do not share memory
        # It is very gross, but I cannot figure out a better way (maybe use install_name_tool to change the dynamic library path to absolute).
        # increment this counter for the hours you wasted on trying find a better way
        # 7 hours

        # this way doesn't work with mulitple isntances fo OVLSolver
        # from . import libavl
        # self.avl = libavl

        module_dir = os.path.dirname(os.path.realpath(__file__))
        module_name = os.path.basename(module_dir)
        if platform.system() == "Windows":
            # HACK
            avl_lib_so_file = glob.glob(os.path.join(module_dir, "libavl*.pyd"))[0]
            print("avl_lib_so_file", avl_lib_so_file)
            if len(avl_lib_so_file) == 0:
                # print the contents of the module dir
                print("module_dir", module_dir)
                print(glob.glob(os.path.join(module_dir, "*")))
                raise RuntimeError("no dyanmic library found")
        else:
            avl_lib_so_file = glob.glob(os.path.join(module_dir, "libavl*.so"))[0]

        # # get just the file name
        avl_lib_so_file = os.path.basename(avl_lib_so_file)

        self.avl = MExt.MExt("libavl", module_name, "optvl", lib_so_file=avl_lib_so_file, debug=debug)._module

        # Initialize AVL
        self.avl.avl()
        if debug:
            self.set_avl_fort_arr("CASE_L", "LVERBOSE", True)

        if timing:
            self.set_avl_fort_arr("CASE_L", "LTIMING", True)

        if geo_file is not None:
            try:
                # check to make sure files exist
                file = geo_file
                f = open(geo_file, "r")
                f.close()

                if mass_file is not None:
                    file = mass_file
                    f = open(mass_file, "r")
                    f.close()
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not open the file '{file}' from python. This is usually an issue with the specified file path"
                )

            self.avl.loadgeo(geo_file)

            if mass_file is not None:
                self.avl.loadmass(mass_file)
        elif input_dict is not None:
            self.load_input_dict(input_dict, postCheck=True)
            self.avl.loadgeo("")
        else:
            raise ValueError("neither a geometry file nor an input options dictionary was specified")

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
            "Cl": 7,
            "Cm": 8,
            "Cl'": 9,
        }

        # control surfaces added in __init__
        # TODO: the keys of this dict aren't used
        self.con_var_to_fort_var = {
            "alpha": ["CASE_R", "ALFA"],
            "beta": ["CASE_R", "BETA"],
        }

        control_names = self.get_control_names()
        self.dindex_to_con_surf = OrderedDict()
        self.con_surf_to_dindex = OrderedDict()
        for idx_c_var, c_name in enumerate(control_names):
            self.dindex_to_con_surf[f"D{idx_c_var + 1}"] = c_name
            self.con_surf_to_dindex[c_name] = f"D{idx_c_var + 1}"

        # set control surface constraint indecies in to con val dict
        idx_control_start = np.max([x for x in self.conval_idx_dict.values()]) + 1
        for idx_c_var, c_name in enumerate(control_names):
            self.conval_idx_dict[c_name] = idx_control_start + idx_c_var
            self.con_var_to_fort_var[c_name] = ["CASE_R", "DELCON"]

        var_to_suffix = {
            "alpha": "AL",
            "beta": "BE",
            "p'": "RX",
            "q'": "RY",
            "r'": "RZ",
        }
        func_to_prefix = {
            "CL": "CL",
            "CD": "CD",
            "CY": "CY",
            "Cl'": "CR",
            "Cm": "CM",
            "Cn'": "CN",
        }
        self.case_stab_derivs_to_fort_var = {
            "spiral parameter": ["CASE_R", "BB"],
            "lateral parameter": ["CASE_R", "RR"],
            "static margin": ["CASE_R", "SM"],
            "neutral point": ["CASE_R", "XNP"],
        }
        for func in func_to_prefix:
            for var in var_to_suffix:
                deriv_key = self._get_deriv_key(var, func)
                self.case_stab_derivs_to_fort_var[deriv_key] = [
                    "CASE_R",
                    f"{func_to_prefix[func]}TOT_{var_to_suffix[var]}",
                ]

        func_to_prefix = {
            "CX": "CX",
            "CY": "CY",
            "CZ": "CZ",
            "Cl": "CR",
            "Cm": "CM",
            "Cn": "CN",
        }

        self.case_body_derivs_to_fort_var = {}
        for func in ["CX", "CY", "CZ", "Cl", "Cm", "Cn"]:
            # for var in var_to_suffix:
            for idx_var, var in enumerate(["u", "v", "w", "p", "q", "r"]):
                deriv_key = self._get_deriv_key(var, func)
                self.case_body_derivs_to_fort_var[deriv_key] = ["CASE_R", f"{func_to_prefix[func]}TOT_U_BA", idx_var]

        #  the case parameters are stored in a 1d array,
        # these indices correspond to the position of each parameter in that arra
        self._init_map_data()

        # set the default solver tolerance
        self.set_avl_fort_arr("CASE_R", "EXEC_TOL", 2e-5)

        if timing:
            print(f"AVL init took {time.time() - start_time} seconds")

    def _init_map_data(self):
        """Used in the __init__ method to allocate the slice data for the surfaces"""
        self.surf_geom_to_fort_var = {}
        self.surf_section_geom_to_fort_var = {}
        self.surf_pannel_to_fort_var = {}
        self.con_surf_to_fort_var = {}
        self.des_var_to_fort_var = {}
        self.body_param_to_fort_var = {}

        self.surface_names = self.get_surface_names()
        self.unique_surface_names = self.get_surface_names(remove_dublicated=True)

        # we have to loop over the unique surfaces because those are the
        # only ones that have geometric data from the input file
        # AVL only mirror the mesh data it doesn't infer the input data
        # for the mirrored surface
        for surf_name in self.unique_surface_names:
            idx_surf = self.get_surface_index(surf_name)

            # only set unduplicated sufaces
            if self.get_avl_fort_arr("SURF_I", "IMAGS", slicer=idx_surf) < 0:
                # this is a duplicated surface, skip it
                raise ValueError("Only non-duplicated surfaces have geom data. Internal list of unique surfaces wrong")

            num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC", slicer=idx_surf)
            nasec_arr = self.get_avl_fort_arr("SURF_GEOM_I", "NASEC", slicer=(idx_surf, slice(None, num_sec)))

            self._setup_surface_maps(surf_name, idx_surf, num_sec)

            self._setup_section_maps(surf_name, idx_surf, num_sec, nasec_arr)

        self.body_names = self.get_body_names()
        self.unique_body_names = self.get_body_names(remove_dublicated=True)

        for body_name in self.unique_body_names:
            idx_body = self.get_body_index(body_name)

            self._setup_body_maps(body_name, idx_body)

    def _setup_surface_maps(self, surf_name: str, idx_surf: int, num_sec: int):
        """Used by the init_map_data and load_input_dict functions to generate which slices of the Fortran array for a
        given geometry, panneling, control surface, or design variable correspond to the given surface. This data is
        stored the surf_geom_to_fort_var dictionary.

        Args:
            surf_name: The name of the surface
            idx_surf: The index of the surface (Fortran indexing)
            num_sec: The number of sections in the surface
        """
        slice_idx_surf = (idx_surf,)
        slice_surf_all = (idx_surf, slice(None))
        slice_surf_secs = (idx_surf, slice(None, num_sec))
        slice_surf_secs_all = slice_surf_secs + (slice(None),)
        slice_surf_secs_x = slice_surf_secs + (0,)
        slice_surf_secs_y = slice_surf_secs + (1,)
        slice_surf_secs_z = slice_surf_secs + (2,)

        self.surf_geom_to_fort_var[surf_name] = {
            "scale": ["SURF_GEOM_R", "XYZSCAL", slice_surf_all],
            "translate": ["SURF_GEOM_R", "XYZTRAN", slice_surf_all],
            "angle": ["SURF_GEOM_R", "ADDINC", slice_idx_surf],
            # "xyzles": ["SURF_GEOM_R", "XYZLES", slice_surf_secs_all],
            "xles": ["SURF_GEOM_R", "XYZLES", slice_surf_secs_x],
            "yles": ["SURF_GEOM_R", "XYZLES", slice_surf_secs_y],
            "zles": ["SURF_GEOM_R", "XYZLES", slice_surf_secs_z],
            "chords": ["SURF_GEOM_R", "CHORDS", slice_surf_secs],
            "aincs": ["SURF_GEOM_R", "AINCS", slice_surf_secs],
            "clcdsec": ["SURF_GEOM_R", "CLCDSEC", slice_surf_secs_all],
            "clcd": ["SURF_GEOM_R", "CLCDSRF", slice_surf_all],  # NEW
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
            "component": ["SURF_I", "LNCOMP", slice_idx_surf],
            "wake": ["SURF_L", "LFWAKE", slice_idx_surf],
            "albe": ["SURF_L", "LFALBE", slice_idx_surf],
            "load": ["SURF_L", "LFLOAD", slice_idx_surf],
        }

        icontd_slices = []
        idestd_slices = []
        xhinged_slices = []
        vhinged_slices = []
        gaind_slices = []
        refld_slices = []
        gaing_slices = []

        for idx_sec in range(num_sec):
            slice_sec = (idx_surf, idx_sec)

            num_con_surf = self.get_avl_fort_arr("SURF_GEOM_I", "NSCON", slicer=slice_sec)
            slice_sec_con = slice_sec + (slice(None, num_con_surf),)
            slice_sec_con_hinge = slice_sec_con + (slice(None),)

            num_des_var = self.get_avl_fort_arr("SURF_GEOM_I", "NSDES", slicer=slice_sec)
            slice_sec_des_var = slice_sec + (slice(None, num_des_var),)

            icontd_slices.append(slice_sec_con)
            xhinged_slices.append(slice_sec_con)
            vhinged_slices.append(slice_sec_con_hinge)
            gaind_slices.append(slice_sec_con)
            refld_slices.append(slice_sec_con)
            idestd_slices.append(slice_sec_des_var)
            gaing_slices.append(slice_sec_des_var)

        self.con_surf_to_fort_var[surf_name] = {
            "icontd": ["SURF_GEOM_I", "ICONTD", icontd_slices],
            "xhinged": ["SURF_GEOM_R", "XHINGED", xhinged_slices],
            "vhinged": ["SURF_GEOM_R", "VHINGED", vhinged_slices],
            "gaind": ["SURF_GEOM_R", "GAIND", gaind_slices],
            "refld": ["SURF_GEOM_R", "REFLD", refld_slices],
        }
        self.des_var_to_fort_var[surf_name] = {
            "idestd": ["SURF_GEOM_I", "IDESTD", idestd_slices],
            "gaing": ["SURF_GEOM_R", "GAING", gaing_slices],
        }

    def _setup_body_maps(self, body_name: str, idx_body: int):
        """Used by the init_map_data and load_input_dict functions to generate which slices of the Fortran array for a
        given geometry or discretization variable correspond to the given body. This data is stored the
        body_param_to_fort_var dictionary.

        Args:
            body_name: The name of the body
            idx_body: The index of the body (Fortran indexing)
        """
        slice_idx_body = (idx_body,)
        slice_body_all = (idx_body, slice(None))

        self.body_param_to_fort_var[body_name] = {
            "nvb": ["BODY_GEOM_I", "NVB", slice_idx_body],
            "bspace": ["BODY_GEOM_R", "BSPACE", slice_idx_body],
            "yduplicate": ["BODY_GEOM_R", "YDUPL_B", slice_idx_body],
            "scale": ["BODY_GEOM_R", "XYZSCAL_B", slice_body_all],
            "translate": ["BODY_GEOM_R", "XYZTRAN_B", slice_body_all],
        }

    def _setup_section_maps(self, surf_name: str, idx_surf: int, num_sec: int, nasec_arr: np.ndarray):
        """Used by the init_map_data and load_input_dict functions to generate which slices of the Fortran array for a
        given section geometry variable correspond to the given surface and section. This data is stored the
        surf_section_geom_to_fort_var dictionary.

        Args:
            surf_name: The name of the surface
            idx_surf: The index of the surface (Fortran indexing)
            num_sec: The number of sections
            nasec_arr: The number of array points used to define the section geometry
        """
        xasec_slices = []
        sasec_slices = []
        casec_slices = []
        tasec_slices = []
        xuasec_slices = []
        xlasec_slices = []
        zuasec_slices = []
        zlasec_slices = []

        for idx_sec in range(num_sec):
            nasec = nasec_arr[idx_sec]
            slice_surf_secs = (idx_surf, idx_sec)
            slice_surf_secs_nasec = slice_surf_secs + (slice(None, nasec),)
            xasec_slices.append(slice_surf_secs_nasec)
            sasec_slices.append(slice_surf_secs_nasec)
            casec_slices.append(slice_surf_secs_nasec)
            tasec_slices.append(slice_surf_secs_nasec)
            xuasec_slices.append(slice_surf_secs_nasec)
            xlasec_slices.append(slice_surf_secs_nasec)
            zuasec_slices.append(slice_surf_secs_nasec)
            zlasec_slices.append(slice_surf_secs_nasec)

        self.surf_section_geom_to_fort_var[surf_name] = {
            # "nasec": [],
            "xasec": ["SURF_GEOM_R", "XASEC", xasec_slices],
            "sasec": ["SURF_GEOM_R", "SASEC", sasec_slices],
            "casec": ["SURF_GEOM_R", "CASEC", casec_slices],
            "tasec": ["SURF_GEOM_R", "TASEC", tasec_slices],
            "xuasec": ["SURF_GEOM_R", "XUASEC", xuasec_slices],
            "xlasec": ["SURF_GEOM_R", "XLASEC", xlasec_slices],
            "zuasec": ["SURF_GEOM_R", "ZUASEC", zuasec_slices],
            "zlasec": ["SURF_GEOM_R", "ZLASEC", zlasec_slices],
        }

    def load_input_dict(self, input_dict: dict, preCheck: bool = True, postCheck: bool = False):
        """Reads and loads the input dictionary data into optvl.
        Equivalent to INPUT routine in AVL.

        Args:
            input_dict: input dictionary in optvl format
            preCheck: perform additional verification of the user's input dictionary before loading into AVL
            postCheck: verify certain inputs values are correctly reflected in the Fortran layer
        """

        # Initialize Variables and Counters
        self.avl.CASE_I.NVOR = 0
        self.avl.CASE_I.NSTRIP = 0
        self.avl.CASE_I.NLNODE = 0
        self.avl.CASE_I.NCONTROL = 0
        self.avl.CASE_I.NDESIGN = 0
        self.set_avl_fort_arr("SURF_GEOM_I", "NSEC", 0, slicer=slice(None, self.NFMAX))
        self.set_avl_fort_arr("BODY_GEOM_I", "NSEC_B", 0, slicer=slice(None, self.NBMAX))
        self.avl.CASE_R.DCL_A0 = 0.0
        self.avl.CASE_R.DCM_A0 = 0.0
        self.avl.CASE_R.DCL_U0 = 0.0
        self.avl.CASE_R.DCM_U0 = 0.0
        self.set_avl_fort_arr("SURF_GEOM_R", "CLCDSRF", 0.0, slicer=(slice(None, 6), slice(None, self.NFMAX)))
        self.avl.CASE_L.LVISC = False
        self.set_avl_fort_arr("SURF_L", "LRANGE", True, slicer=slice(None, self.NFMAX))

        # Perform pre-check of user's input dictionary before loading into AVL
        if preCheck:
            input_dict = pre_check_input_dict(input_dict)

        def get_types_from_blk(common_blk):
            """Determines the variable type a common block uses from
            its name

            Args:
                common_blk: Name of the common block

            Returns:
                Python type equivalent used by the common block
            """
            last_char = common_blk[-1]

            if last_char == "C":
                return str
            elif last_char == "R":
                return (float, np.float64)
            elif last_char == "I":
                return (int, np.int32)
            elif last_char == "L":
                return (bool, int, np.int32)
            else:
                raise ValueError(f"type not able to be infered from common block {common_blk}")

        def check_type(key, avl_vars, given_val):
            """Checks the type for a given AVL Fortran Common Block var against a given value

            Args:
                key: OptVL input variable dictionary key
                avl_vars: AVL Common Block variable name array
                given_val: Input value to check type against

            """
            # get the type that it should be
            expected_type = get_types_from_blk(avl_vars[0])

            # if the expected type is a str
            if expected_type is str:
                # check the type of the scaler
                if not isinstance(given_val, expected_type):
                    raise TypeError(f"Variable {key} is of type {type(given_val)} but expected {expected_type}")

                # for strings no further checks are required
                return

            # --- test shape and type of numeric vales---
            current_val = self.get_avl_fort_arr(*avl_vars)

            # is the current value a scalar or a numpy array?
            if isinstance(current_val, np.ndarray):
                if not isinstance(given_val, np.ndarray):
                    raise ValueError(
                        f"Variable {key} is scalar, but optvl expected an array of shape {current_val.shape}"
                    )

                # compare the shapes
                if current_val.shape != val.shape:
                    raise ValueError(
                        f"Variable {key} is shape {given_val.shape}, but optvl expected an array of shape {current_val.shape}"
                    )

                # check that the type of the array matches the expectation
                if not isinstance(given_val.flatten()[0], expected_type):
                    raise TypeError(
                        f"Variable {key} is an array of type {given_val.dtype} but expected {expected_type}"
                    )
            else:
                # check the type of the scaler
                if not isinstance(given_val, expected_type):
                    raise TypeError(
                        f"Variable {key} is a scalar of type {type(given_val)} but expected {expected_type}"
                    )

        # Set AVL header variables
        # CDp is the only optional input for the AVL header
        optional_header_defaults = {"CDp": 0.0}

        for key, avl_vars in self.header_var_to_fort_var.items():
            if key not in input_dict:
                if key in optional_header_defaults:
                    val = optional_header_defaults[key]
                else:
                    raise ValueError(f"Key {key} not found in input dictionary but is required")
            else:
                val = input_dict[key]

            check_type(key, avl_vars, val)

            self.set_avl_fort_arr(avl_vars[0], avl_vars[1], val)

        self.set_avl_fort_arr("CASE_R", "YSYM", 0.0)  # YSYM Hardcoded to 0

        # set the global control variable options
        ncontrol = len(input_dict.get("dname", []))
        if ncontrol > self.NDMAX:
            raise RuntimeError(f"Number of specified controls exceeds {self.NDMAX}. Raise NDMAX!")
        self.set_avl_fort_arr("CASE_I", "NCONTROL", ncontrol)

        for k in range(ncontrol):
            self.avl.CASE_C.DNAME[k] = input_dict["dname"][k]

        # set the gloabl design variable options
        ndesign = len(input_dict.get("gname", []))
        self.set_avl_fort_arr("CASE_I", "NDESIGN", ncontrol)
        if ndesign > self.NGMAX:
            raise RuntimeError(f"Number of specified design variables exceeds {self.NGMAX}. Raise NGMAX!")

        for k in range(ncontrol):
            self.avl.CASE_C.GNAME[k] = input_dict["gname"][k]

        # Set total number of surfaces in one shot
        num_surfs = len(input_dict["surfaces"])
        if num_surfs < self.NFMAX:
            self.set_avl_fort_arr("CASE_I", "NSURF", num_surfs)  # YSYM Hardcoded to 0
        else:
            raise RuntimeError(f"Number of specified surfaces, {num_surfs}, exceeds {self.NFMAX}. Raise NFMAX!")

        # Load surfaces
        if num_surfs > 0:
            surf_names = list(input_dict["surfaces"].keys())

            # setup surface data for initial input
            self.surf_geom_to_fort_var = {}
            self.surf_section_geom_to_fort_var = {}
            self.surf_pannel_to_fort_var = {}
            self.con_surf_to_fort_var = {}
            self.des_var_to_fort_var = {}

            idx_surf = 0

            for surf_name in input_dict["surfaces"]:
                surf_dict = input_dict["surfaces"][surf_name]

                num_secs = surf_dict["num_sections"]
                # Set total number of sections in one shot
                if num_secs < self.NSMAX:
                    self.set_avl_fort_arr("SURF_GEOM_I", "NSEC", num_secs, slicer=idx_surf)
                else:
                    raise RuntimeError(
                        f"Number of specified sections for surface {surf_name} exceeds {self.NSMAX}. Raise NSMAX!"
                    )

                # Set the number of control and design variables for the surface
                for idx_sec in range(num_secs):
                    self.set_avl_fort_arr(
                        "SURF_GEOM_I", "NSCON", surf_dict["num_controls"][idx_sec], slicer=(idx_surf, idx_sec)
                    )
                    self.set_avl_fort_arr(
                        "SURF_GEOM_I", "NSDES", surf_dict["num_design_vars"][idx_sec], slicer=(idx_surf, idx_sec)
                    )

                self._setup_surface_maps(surf_name, idx_surf, num_secs)

                # Set surface name
                self.avl.CASE_C.STITLE[idx_surf] = self._str_to_fort_str(surf_name, num_max_char=40)

                # fmt: off

                optional_surface_defaults = {
                    "nspan":                  0,                      
                    "sspace":                 0.0,                    
                    "use surface spacing":   False,                   
                    "component":              idx_surf+1, # +1 for 1-based indexing in fortran
                    "scale":                  np.array([1.,1.,1.]),   
                    "translate":              np.array([0.,0.,0.]),   
                    "angle":                  0.0,                    
                    "wake":                   True,                   
                    "albe":                   True,                   
                    "load":                   True,                   
                    "clcd":                   np.zeros(6, dtype=np.float64),
                    "nspans":                 np.zeros(num_secs, dtype=int),
                    "sspaces":                np.zeros(num_secs, dtype=int),
                    "clcdsec":                np.zeros((num_secs,6)),
                    "claf":                   np.ones(num_secs),
                }
                # fmt: on

                # set some flags based on the options used for this surface
                if "sspace" in surf_dict:
                    self.set_avl_fort_arr("SURF_GEOM_L", "LSURFSPACING", True, slicer=idx_surf)
                else:
                    self.set_avl_fort_arr("SURF_GEOM_L", "LSURFSPACING", False, slicer=idx_surf)

                if "yduplicate" in surf_dict:
                    self.set_avl_fort_arr("SURF_GEOM_L", "LDUPL", True, slicer=idx_surf)
                else:
                    self.set_avl_fort_arr("SURF_GEOM_L", "LDUPL", False, slicer=idx_surf)

                if "clcd" in surf_dict or "clcdsec" in surf_dict:
                    # if any of the surfaces use clcd then turn on viscous loads
                    self.set_avl_fort_arr("CASE_L", "LVISC", True)

                # lhinge = False Appears in AVL but does nothing

                for key, avl_vars in chain(
                    self.surf_geom_to_fort_var[surf_name].items(), self.surf_pannel_to_fort_var[surf_name].items()
                ):
                    if key not in surf_dict:
                        if key in optional_surface_defaults:
                            val = optional_surface_defaults[key]
                        else:
                            raise ValueError(f"Key {key} not found in surface dictionary, {surf_name}, but is required")
                    else:
                        val = surf_dict[key]

                    check_type(key, avl_vars, val)

                    self.set_avl_fort_arr(avl_vars[0], avl_vars[1], val, slicer=avl_vars[2])

                # determine what method of airfoil definition we are using
                # check to make sure we don't have multiple airfoil definitions used for this surface
                airfoil_spec_keys = {"naca", "airfoils", "afiles", "xasec"} & surf_dict.keys()

                if len(airfoil_spec_keys) > 1:
                    raise KeyError(
                        f"OptVL can only have one method of specifing airfoil geometry per surface, found {airfoil_spec_keys} in surface {surf_name}"
                    )

                xfminmax_arr = surf_dict.get("xfminmax", np.array([0.0, 1.0] * num_secs))
                num_pts = min(50, self.IBX)

                # setup for manually specifying coordinates
                if "xasec" in surf_dict.keys():
                    warnings.warn(
                        "OptVL WARNING - Setting airfoil section data directly via the input dictionary is not recommened!\n "
                        "OptVL will not verify that the camber line is consistent with the given coordiantes.\n"
                        "Specify the coordinates directly in the inputs dictionary with ''airfoils'' or ''afiles'' or use the set_section_coordinates function.",
                        stacklevel=2,
                    )

                    # do we have section data?
                    nasec_list = [len(x) for x in surf_dict["xasec"]]
                    self._setup_section_maps(surf_name, idx_surf, num_secs, nasec_list)

                # Load airfoil data sections
                # Load the Airfoil Section into AVL
                for j in range(num_secs):
                    xfminmax = xfminmax_arr[j]

                    # Manually Specify Coordiantes (no camberline verification, only use if you know what you're doing)
                    if "xasec" in surf_dict.keys():
                        self.set_avl_fort_arr("SURF_GEOM_I", "NASEC", nasec_list[j], slicer=(idx_surf, j))

                        # TODO-JLA: a user should not have to spcify XLASEC, XUASEC, ZLASEC, ZUASEC
                        # since these can be calculated from the other inputs

                        for key in self.surf_section_geom_to_fort_var[surf_name]:
                            avl_vars_secs = self.surf_section_geom_to_fort_var[surf_name][key]
                            avl_vars = (avl_vars_secs[0], avl_vars_secs[1], avl_vars_secs[2][j])

                            if key not in surf_dict:
                                raise ValueError(
                                    f"Key `{key}` not found in surface dictionary, `{surf_name}`, but is required when manually specifing airfoil coordinates"
                                )

                            val = surf_dict[key][j]

                            check_type(key, avl_vars, val)
                            self.set_avl_fort_arr(avl_vars[0], avl_vars[1], val, slicer=avl_vars[2])

                    # 4 digit NACA airfoil specification
                    elif "naca" in surf_dict.keys():
                        if (xfminmax[0] > 0.01) or (xfminmax[1] < 0.99):
                            self.set_avl_fort_arr("SURF_L", "LRANGE", True, slicer=idx_surf)

                        # Store this stuff so we can read it later
                        self.avl.CASE_C.NACA[j, idx_surf] = surf_dict["naca"][j]
                        self.avl.SURF_GEOM_R.XFMIN_R[j, idx_surf] = xfminmax[0]
                        self.avl.SURF_GEOM_R.XFMAX_R[j, idx_surf] = xfminmax[1]
                        self.set_section_naca(j, idx_surf, num_pts, surf_dict["naca"][j], xfminmax)

                    # Airfoil coordinates set directly in dictionary
                    elif "airfoils" in surf_dict.keys():
                        self.set_section_coordinates(
                            j,
                            idx_surf,
                            num_pts,
                            surf_dict["airfoils"][j][0],
                            surf_dict["airfoils"][j][1],
                            xfminmax,
                            storecoords=True,
                        )

                    # Load airfoil file
                    elif "afiles" in surf_dict.keys():
                        self.avl.CASE_C.AFILES[j, idx_surf] = surf_dict["afiles"][j]
                        X = read_coordinates_files(surf_dict["afiles"][j])
                        self.set_section_coordinates(j, idx_surf, num_pts, X[:, 0], X[:, 1], xfminmax, storecoords=True)

                    else:
                        # apply a flat camber line
                        slicer_airfoil_flat = (idx_surf, j, slice(0, 2))
                        self.set_avl_fort_arr("SURF_GEOM_I", "NASEC", 2, slicer=(idx_surf, j))
                        self.set_avl_fort_arr("SURF_GEOM_R", "XASEC", np.array([0.0, 1.0]), slicer=slicer_airfoil_flat)
                        self.set_avl_fort_arr("SURF_GEOM_R", "SASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)
                        self.set_avl_fort_arr("SURF_GEOM_R", "TASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)

                        self.set_avl_fort_arr("SURF_GEOM_R", "XLASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)
                        self.set_avl_fort_arr("SURF_GEOM_R", "XUASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)
                        self.set_avl_fort_arr("SURF_GEOM_R", "ZLASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)
                        self.set_avl_fort_arr("SURF_GEOM_R", "ZUASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)
                        self.set_avl_fort_arr("SURF_GEOM_R", "CASEC", np.array([0.0, 0.0]), slicer=slicer_airfoil_flat)

                # --- setup control variables for each section ---
                # Load control surfaces
                if "icontd" in surf_dict.keys():
                    for j in range(num_secs):
                        # check to make sure this section has control vars
                        if surf_dict["num_controls"][j] == 0:
                            continue

                        for key, avl_vars in self.con_surf_to_fort_var[surf_name].items():
                            avl_vars_secs = self.con_surf_to_fort_var[surf_name][key]
                            avl_vars = (avl_vars_secs[0], avl_vars_secs[1], avl_vars_secs[2][j])

                            if key not in surf_dict:
                                raise ValueError(
                                    f"Key {key} not found in surf dictionary, `{surf_name}` but is required"
                                )
                            else:
                                val = surf_dict[key][j]

                            check_type(key, avl_vars, val)
                            self.set_avl_fort_arr(avl_vars[0], avl_vars[1], val, slicer=avl_vars[2])

                # --- setup design variables for each section ---
                # Load design variables
                if "idestd" in surf_dict.keys():
                    for j in range(num_secs):
                        # check to make sure this section has control vars
                        if surf_dict["num_design_vars"][j] == 0:
                            continue

                        for key, avl_vars in self.des_var_to_fort_var[surf_name].items():
                            avl_vars_secs = self.des_var_to_fort_var[surf_name][key]
                            avl_vars = (avl_vars_secs[0], avl_vars_secs[1], avl_vars_secs[2][j])

                            if key not in surf_dict:
                                raise ValueError(
                                    f"Key {key} not found in surf dictionary, `{surf_name}` but is required"
                                )
                            else:
                                val = surf_dict[key][j]

                            check_type(key, avl_vars, val)
                            self.set_avl_fort_arr(avl_vars[0], avl_vars[1], val, slicer=avl_vars[2])

                # Make the surface
                if self.debug:
                    print(f"Building surface: {surf_name}")
                self.avl.makesurf(idx_surf + 1)  # +1 to convert to 1 based indexing

                if "yduplicate" in surf_dict.keys():
                    self.avl.sdupl(idx_surf + 1, surf_dict["yduplicate"], "YDUP")
                    self.avl.CASE_I.NSURF += 1
                    idx_surf += 1

                    # Keep python data consistent with Fortran
                    surf_names.insert(idx_surf + 1, surf_name)

                idx_surf += 1

        # Set total number of bodies in one shot
        if len(input_dict["bodies"]) < self.NBMAX:
            self.avl.CASE_I.NBODY = len(input_dict["bodies"])
        else:
            raise RuntimeError(f"Number of specified surfaces exceeds {self.NBMAX}. Raise NBMAX!")

        # Load bodies
        if len(input_dict["bodies"]) > 0:
            self.body_param_to_fort_var = {}

            idx_body = 0

            for body_name in input_dict["bodies"]:
                body_dict = input_dict["bodies"][body_name]

                self._setup_body_maps(body_name, idx_body)

                # Setup body
                self.avl.BODY_GEOM_I.NSEC_B[idx_body] = 0  # Body sections not used but set to 0 just in case

                # Set body name
                self.avl.CASE_C.BTITLE[idx_body] = self._str_to_fort_str(body_name, num_max_char=40)

                if "yduplicate" in body_dict:
                    self.set_avl_fort_arr("BODY_GEOM_L", "LDUPL_B", True, slicer=idx_body)
                else:
                    self.set_avl_fort_arr("BODY_GEOM_L", "LDUPL_B", False, slicer=idx_body)

                optional_body_defaults = {"yduplicate": 0.0}

                for key, avl_vars in self.body_param_to_fort_var[body_name].items():
                    if key not in body_dict:
                        if key in optional_body_defaults:
                            val = optional_body_defaults[key]
                        else:
                            raise ValueError(f"Key {key} not found in body dictionary, {body_name}, but is required")
                    else:
                        val = body_dict[key]

                    check_type(key, avl_vars, val)

                    self.set_avl_fort_arr(avl_vars[0], avl_vars[1], val, slicer=avl_vars[2])

                # Load airfoil file
                if "body_oml" in body_dict.keys():
                    self.set_body_coordinates(
                        idx_body, num_pts, body_dict["body_oml"][0, :], body_dict["body_oml"][1, :], True
                    )
                if "bfile" in body_dict.keys():
                    # xfminmax doesn't appear to be supported for bodies and the fact that the bfil input routine in the fortran layer takes them at all is a bug?
                    # xfminmax = body_dict["xfminmax"][j] if "xfminmax" in body_dict.keys() else np.array([0., 1.])
                    # if ((xfminmax[0] > 0.01) or (xfminmax[1] < 0.99)):
                    #     self.set_avl_fort_arr("SURF_L","LRANGE", True, slicer=i)
                    self.avl.CASE_C.BFILES[idx_body] = body_dict["bfile"]
                    X = read_coordinates_files(body_dict["bfile"])
                    self.set_body_coordinates(idx_body, num_pts, X[:, 0], X[:, 1], True)

                # Make the body
                if self.debug:
                    print(f"Building body: {body_name}")
                self.avl.makebody(idx_body + 1)

                if "yduplicate" in body_dict.keys():
                    self.avl.bdupl(idx_body + 1, body_dict["yduplicate"], "YDUP")

                    # Keep python data consistent with Fortran
                    idx_body += 1

                idx_body += 1

        if postCheck:
            self.post_check_input(input_dict)

        if self.debug:
            print(f"Mach: {self.avl.CASE_R.MACH0}")
            print(f"Number of Bodies: {self.avl.CASE_I.NBODY}")
            print(f"Number of Surfaces: {self.avl.CASE_I.NSURF}")
            print(f"Number of Strips: {self.avl.CASE_I.NSTRIP}")
            print(f"Number of Horseshoe Vortices: {self.avl.CASE_I.NVOR}")
            print(f"Number of Control Surfaces: {self.avl.CASE_I.NCONTROL}")
            print(f"Number of Design Parameters (AVL): {self.avl.CASE_I.NDESIGN}")

            if self.avl.CASE_I.IYSYM > 0:
                print(f"Y Symmetry: Wall plane  at Ysym = {self.avl.CASE_R.YSYM}")
            elif self.avl.CASE_I.IYSYM < 0:
                print(f"Y Symmetry: Free surface at Ysym = {self.avl.CASE_R.YSYM}")

            if self.avl.CASE_I.IZSYM > 0:
                print(f"Z Symmetry: Ground plane at Zsym = {self.avl.CASE_R.ZSYM}")
            elif self.avl.CASE_I.IZSYM < 0:
                print(f"Z Symmetry: Free surface at Zsym = {self.avl.CASE_R.ZSYM}")

        # Tell AVL that geometry exists now and is ready for analysis
        self.avl.CASE_L.LGEO = True

    def set_section_naca(self, isec: int, isurf: int, nasec: int, naca: str, xfminmax: np.ndarray):
        """Sets the airfoil oml points for the specified surface and section. Computes camber lines, thickness, and oml shape from
        NACA 4-digit specification.


        Args:
            isec: section number to set the airfoil mesh
            isurf: surface number to set the airfoil mesh
            nasec: number of points to evaluate the interpolated camber line and thickness curves at
            naca: 4-digit naca specificaion as a string
            xfminmax: length 2 array with the min and max x/c to slice the airfoil
        """

        # Read 4-digit string
        cam = np.float64(naca[0]) / 100.0
        pos = np.float64(naca[1]) / 10.0
        thick = np.float64(naca[2:]) / 100.0

        # Generate airfoil section data
        xf = xfminmax[0] + np.diff(xfminmax) * np.arange(nasec) / (nasec - 1)
        slopes = np.zeros_like(xf)
        slopes[xf < pos] = 2.0 * cam * (pos - xf[xf < pos]) / (pos**2)
        slopes[xf > pos] = 2.0 * cam * (pos - xf[xf > pos]) / (1 - pos) ** 2
        thickness = (
            10.0 * thick * (0.29690 * np.sqrt(xf) - 0.12600 * xf - 0.35160 * xf**2 + 0.28430 * xf**3 - 0.10150 * xf**4)
        )
        zf = np.zeros_like(xf)
        zf[xf < pos] = cam * (2.0 * pos * xf[xf < pos] - 1.0) * xf[xf < pos] / (pos**2)
        zf[xf > pos] = cam * ((1 - 2.0 * pos) + (2.0 * pos - xf[xf > pos]) * xf[xf > pos]) / (1 - pos) ** 2
        theta = np.arctan(slopes)

        # Set airfoil section
        self.set_avl_fort_arr("SURF_GEOM_I", "NASEC", nasec, slicer=(isurf, isec))
        self.set_avl_fort_arr(
            "SURF_GEOM_R", "XASEC", (xf - xf[0]) / (xf[-1] - xf[0]), slicer=(isurf, isec, slice(0, nasec))
        )
        self.set_avl_fort_arr("SURF_GEOM_R", "SASEC", slopes, slicer=(isurf, isec, slice(0, nasec)))
        self.set_avl_fort_arr("SURF_GEOM_R", "TASEC", thickness, slicer=(isurf, isec, slice(0, nasec)))

        self.set_avl_fort_arr(
            "SURF_GEOM_R", "XLASEC", xf + 0.5 * thickness * np.sin(theta), slicer=(isurf, isec, slice(0, nasec))
        )
        self.set_avl_fort_arr(
            "SURF_GEOM_R", "XUASEC", xf - 0.5 * thickness * np.sin(theta), slicer=(isurf, isec, slice(0, nasec))
        )
        self.set_avl_fort_arr(
            "SURF_GEOM_R", "ZLASEC", zf - 0.5 * thickness * np.cos(theta), slicer=(isurf, isec, slice(0, nasec))
        )
        self.set_avl_fort_arr(
            "SURF_GEOM_R", "ZUASEC", zf + 0.5 * thickness * np.cos(theta), slicer=(isurf, isec, slice(0, nasec))
        )
        self.set_avl_fort_arr("SURF_GEOM_R", "CASEC", zf, slicer=(isurf, isec, slice(0, nasec)))

    def set_section_coordinates(
        self,
        isec: int,
        isurf: int,
        nasec: int,
        x: np.ndarray,
        y: np.ndarray,
        xfminmax: np.ndarray,
        storecoords: bool = False,
    ):
        """Sets the airfoil oml points for the specified surface and section. Computes the camber line and interpolates it
        with AVL's 1D Akima Spline implementation.


        Args:
            isec: section number to set the airfoil mesh
            isurf: surface number to set the airfoil mesh
            nasec: number of points to evaluate the interpolated camber line and thickness curves at
            x: airfoil x-coordinate array
            y: airfoil y-coodinate array
            xfminmax: length 2 array with the min and max x/c to slice the airfoil
            storecoords: store the raw input coordinates in common block
        """

        if (len(x) > self.IBX) or (len(y) > self.IBX):
            raise RuntimeError(f"Airfoil array overflow! Increase IBX to {len(x) if len(x) > len(y) else len(y)}")
        if (len(x) < 2) or (len(y) < 2):
            raise RuntimeError("Airfoil shape not defined. Too few points!")
        if len(x) != len(y):
            raise RuntimeError(f"x and y array lengths do not match! len x: {len(x)}, len y: {len(y)}")
        if isurf + 1 > self.get_num_surfaces():
            raise RuntimeError(f"surface {isurf} does not exist!")
        # if isec+1 > self.get_num_sections(self.get_surface_names(remove_dublicated=True)[isec]):
        #     raise RuntimeError(f"section {isec} in surface {isurf} does not exist!")

        self.avl.set_section_coordinates(isec + 1, isurf + 1, x, y, nasec, xfminmax[0], xfminmax[1], storecoords)

    def set_body_coordinates(self, ibod: int, nasec: int, x: np.ndarray, y: np.ndarray, storecoords: bool = False):
        """Sets the body of revolution oml points for the specified body. Computes the camber line and interpolates it
        with AVL's 1D Akima Spline implementation.


        Args:
            ibod: body number to set the outer mold line too
            nasec: number of points to evaluate the interpolated camber line and thickness curves at
            x: oml x-coordinate array
            y: oml y-coodinate array
            xfminmax: length 2 array with the min and max x/c to slice the oml
            storecoords: store the raw input coordinates in common block
        """

        if (len(x) > self.IBX) or (len(y) > self.IBX):
            raise RuntimeError(f"Body oml array overflow! Increase IBX to {len(x) if len(x) > len(y) else len(y)}")
        if (len(x) < 2) or (len(y) < 2):
            raise RuntimeError("Airfoil shape not defined. Too few points!")
        if len(x) != len(y):
            raise RuntimeError(f"x and y array lengths do not match! len x: {len(x)}, len y: {len(y)}")
        if ibod + 1 > self.avl.CASE_I.NBODY:
            raise RuntimeError(f"body {ibod} does not exist!")
        # if isec+1 > self.get_num_sections(self.get_surface_names(remove_dublicated=True)[isec]):
        #     raise RuntimeError(f"section {isec} in surface {isurf} does not exist!")

        self.avl.set_body_coordinates(ibod + 1, x, y, nasec, storecoords)

    # MOVED TO FORTRAN...amake.f
    # def set_section_coordinates(self,isec: int, isurf: int, nasec: int, x: np.ndarray, y: np.ndarray, xfminmax: np.ndarray):
    #     """Sets the airfoil oml points for the specified surface and section. Computes the camber line and interpolates it
    #     with AVL's 1D Akima Spline implementation.

    #     Args:
    #         isec: section number to set the airfoil mesh
    #         isurf: surface number to set the airfoil mesh
    #         nasec: number of points to evaluate the interpolated camber line and thickness curves at
    #         x: airfoil x-coordinate array
    #         y: airfoil y-coodinate array
    #         xfminmax: length 2 array with the min and max x/c to slice the airfoil
    #     """

    #     if ((len(x) > self.IBX) or (len(y) > self.IBX)):
    #         raise RuntimeError(f"Airfoil array overflow! Increase IBX to {len(x) if len(x)>len(y) else len(y)}")
    #     if len(x) != len(y):
    #         raise RuntimeError(f"x and y array lengths do not match! len x: {len(x)}, len y: {len(y)}")

    #     xf = xfminmax[0] + np.diff(xfminmax)*np.arange(nasec)/(nasec-1)

    #     xin,yin,tin,nasec = self.avl.getcam(x,y,len(x),nasec,True)

    #     xasec = xin[0] + xf*(xin[-1] -xin[0])

    #     zc = np.zeros(nasec)
    #     slopes = np.zeros(nasec)
    #     thickness = np.zeros(nasec)

    #     # Use the original AVL akima, not vectorized
    #     for i in range(nasec):
    #         zc[i], slopes[i] = self.avl.akima(xin, yin, nasec, xasec) # USE_CPOML always on in OptVL
    #         thickness[i], _ = self.avl.akima(xin, tin, nasec, xasec)

    #     # Set airfoil section
    #     self.set_avl_fort_arr("SURF_GEOM_I","NASEC", nasec, slicer=(isurf,isec))
    #     self.set_avl_fort_arr("SURF_GEOM_R","XASEC", (xasec-xasec[0])/(xasec[-1]-xasec[0]), slicer=(isurf,isec,slice(0,nasec)))
    #     self.set_avl_fort_arr("SURF_GEOM_R","SASEC", slopes, slicer=(isurf,isec,slice(0,nasec)))
    #     self.set_avl_fort_arr("SURF_GEOM_R","TASEC", thickness, slicer=(isurf,isec,slice(0,nasec)))

    #     self.set_avl_fort_arr("SURF_GEOM_R","XLASEC", xasec, slicer=(isurf,isec,slice(0,nasec)))
    #     self.set_avl_fort_arr("SURF_GEOM_R","XUASEC", xasec, slicer=(isurf,isec,slice(0,nasec)))
    #     self.set_avl_fort_arr("SURF_GEOM_R","ZLASEC", zc - 0.5*thickness, slicer=(isurf,isec,slice(0,nasec)))
    #     self.set_avl_fort_arr("SURF_GEOM_R","ZUASEC", zc + 0.5*thickness, slicer=(isurf,isec,slice(0,nasec)))
    #     self.set_avl_fort_arr("SURF_GEOM_R","CASEC", zc, slicer=(isurf,isec,slice(0,nasec)))

    def post_check_input(self, inputDict: dict):
        """This routine verifies that a few critical values in the Fortran layer have been set correctly with
        regard to the input dict.

        To be expanded later...
        """

        # Surface checks
        if "surfaces" in inputDict.keys():
            # check number of surfaces
            dict_num_surfaces = len(inputDict["surfaces"])
            surf_names = list(inputDict["surfaces"].keys())
            for i in range(len(inputDict["surfaces"])):
                surf_dict = inputDict["surfaces"][surf_names[i]]
                if "yduplicate" in surf_dict.keys():
                    dict_num_surfaces += 1
            if dict_num_surfaces != self.avl.CASE_I.NSURF:
                raise RuntimeError(
                    f"Mismatch: NSURF = {self.avl.CASE_I.NSURF}, Dictionary: {len(inputDict['surfaces'])}"
                )

            # check number of sections, controls, and dvs
            if len(inputDict["surfaces"]) > 0:
                surf_names = list(inputDict["surfaces"].keys())
                for i in range(len(inputDict["surfaces"])):
                    surf_dict = inputDict["surfaces"][surf_names[i]]
                    if self.avl.SURF_GEOM_I.NSEC[i] != surf_dict["num_sections"]:
                        raise RuntimeError(
                            f"Mismatch: NSEC[i] = {self.avl.SURF_GEOM_I.NSEC[i]}, Dictionary: {surf_dict['num_sections']}"
                        )

                    # Check controls and design variables per section
                    for j in range(surf_dict["num_sections"]):
                        if self.avl.SURF_GEOM_I.NSCON[i, j] != surf_dict["num_controls"][j]:
                            raise RuntimeError(
                                f"Mismatch: NSCON[i,j] = {self.avl.SURF_GEOM_I.NSCON[i, j]}, Dictionary: {surf_dict['num_controls'][j]}"
                            )
                        if self.avl.SURF_GEOM_I.NSDES[i, j] != surf_dict["num_design_vars"][j]:
                            raise RuntimeError(
                                f"Mismatch: NSDES[i,j] = {self.avl.SURF_GEOM_I.NSDES[i, j]}, Dictionary: {surf_dict['num_design_vars'][j]}"
                            )

                # Check the global control and design var count
                if "dname" in inputDict.keys():
                    if len(inputDict["dname"]) != self.avl.CASE_I.NCONTROL:
                        raise RuntimeError(
                            f"Mismatch: NCONTROL = {self.avl.CASE_I.NCONTROL}, Dictionary: {inputDict['dname']}"
                        )
                if "gname" in inputDict.keys():
                    if len(inputDict["gname"]) != self.avl.CASE_I.NDESIGN:
                        raise RuntimeError(
                            f"Mismatch: NDESIGN = {self.avl.CASE_I.NDESIGN}, Dictionary: {inputDict['gname']}"
                        )

        # Body checks
        if "bodies" in inputDict.keys():
            dict_num_bodies = len(inputDict["bodies"])
            body_names = list(inputDict["bodies"].keys())
            for i in range(len(inputDict["bodies"])):
                body = inputDict["bodies"][body_names[i]]
                if "yduplicate" in body.keys():
                    dict_num_bodies += 1
            # check number of bodies
            if dict_num_bodies != self.avl.CASE_I.NBODY:
                raise RuntimeError(f"Mismatch: NBODY = {self.avl.CASE_I.NBODY}, Dictionary: {len(inputDict['bodies'])}")

    # region -- analysis api
    def execute_run(self, tol: float = 0.00002):
        """Run the analysis (equivalent to the AVL command `x` in the OPER menu)

        Args:
            tol: the tolerace of the Newton solver used for timing the aircraft
        """
        self.set_avl_fort_arr("CASE_R", "EXEC_TOL", tol)
        self.avl.oper()

    def set_variable(self, var: str, val: float):
        """set a variable for the run case (equivalent to setting a variable in AVL's OPER menu)
        Args:
            var: variable to be constrained ["alpha"", "beta"", "roll rate", "pitch rate", "yaw rate"] or any control surface.
            val: target value of `con_var`
        """
        avl_variables = {
            "alpha": ("CASE_R", "ALFA"),
            "beta": ("CASE_R", "BETA"),
            "roll rate": ("CASE_R", "WROT", 0),
            "pitch rate": ("CASE_R", "WROT", 1),
            "yaw rate": ("CASE_R", "WROT", 2),
        }
        ref_data = self.get_reference_data()
        dtr = self.get_avl_fort_arr("CASE_R", "DTR")
        bref = ref_data["Bref"]
        cref = ref_data["Cref"]

        vars_factors = {
            "alpha": dtr,
            "beta": dtr,
            "roll rate": 2 / bref,
            "pitch rate": 2 / cref,
            "yaw rate": 2 / bref,
        }

        if var in avl_variables:
            # save the name of the avl_var 0
            avl_var = avl_variables[var]
        else:
            raise ValueError(
                f"specified variable `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.con_surf_to_dindex.items()]}. Constraints that must be implicitly satisfied (such as `CL`) are set with `add_trim_constraint`."
            )

        # check that the type of val is correct
        if isinstance(val, (int, float, np.floating, np.integer)):
            pass
        elif isinstance(val, np.ndarray):
            if val.size != 1:
                raise TypeError(f"variable value must be can only be an np.ndarray if it has size = 1. Got {val}")
        else:
            raise TypeError(f"variable value must be a int or float for contraint {var}. Got {val}")

        slicer = (avl_var[2],) if len(avl_var) == 3 else None
        self.set_avl_fort_arr(avl_var[0], avl_var[1], val * vars_factors[var], slicer)

        # set the contraint value so that the set value is used in analysis
        self.set_avl_fort_arr("CASE_R", "CONVAL", val, (self.conval_idx_dict[var],))

    def get_variable(self, var: str):
        """set a variable for the run case (equivalent to setting a variable in AVL's OPER menu)
        Args:
            var: variable to be constrained ["alpha"", "beta"", "roll rate", "pitch rate", "yaw rate"] or any control surface.
            val: target value of `con_var`
        """
        avl_variables = {
            "alpha": ("CASE_R", "ALFA"),
            "beta": ("CASE_R", "BETA"),
            "roll rate": ("CASE_R", "WROT", 0),
            "pitch rate": ("CASE_R", "WROT", 1),
            "yaw rate": ("CASE_R", "WROT", 2),
        }

        ref_data = self.get_reference_data()
        dtr = self.get_avl_fort_arr("CASE_R", "DTR")
        bref = ref_data["Bref"]
        cref = ref_data["Cref"]

        vars_factors = {
            "alpha": dtr,
            "beta": dtr,
            "roll rate": 2 / bref,
            "pitch rate": 2 / cref,
            "yaw rate": 2 / bref,
        }

        # Check that the variable is valid
        if var in avl_variables:
            # save the name of the avl_var
            var_key = avl_variables[var]
        else:
            raise ValueError(
                f"specified variable `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.con_surf_to_dindex.items()]}."
            )

        slicer = var_key[2] if len(var_key) == 3 else None

        val = copy.deepcopy(self.get_avl_fort_arr(var_key[0], var_key[1], slicer))
        val /= vars_factors[var]

        return val

    def set_constraint(self, var: str, con_var: str, val: float):
        """Set the constraints on the analysis case (equivalent to setting a constraint in AVL's OPER menu)

        Args:
            var: variable to be constrained ["alpha"", "beta"", "roll rate", "pitch rate", "yaw rate"] or any control surface.
            val: target value of `con_var`
            con_var: variable output that needs to be constrained. It could be any value for `var` plus ["CL", "CY", "Cl", "Cm", "Cn"]. If None, than `var` is also the `con_var`

        """
        avl_variables = {
            "alpha": "A ",
            "beta": "B ",
            "roll rate": "R ",
            "pitch rate": "P ",
            "yaw rate": "Y ",
        }

        avl_con_variables = copy.deepcopy(avl_variables)
        avl_con_variables.update(
            {
                "CL": "C ",
                "CY": "S ",
                "Cl": "RM",
                "Cm": "PM",
                "Cn": "YM",
            }
        )

        if var in avl_variables:
            # save the name of the avl_var
            avl_var = avl_variables[var]
        elif var in self.con_surf_to_dindex:
            avl_var = self.con_surf_to_dindex[var]
        elif var in self.dindex_to_con_surf:
            avl_var = var
        else:
            raise ValueError(
                f"specified key `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.con_surf_to_dindex.items()]}."
            )

        if con_var is None:
            avl_con_var = avl_var
        elif con_var in avl_con_variables:
            avl_con_var = avl_con_variables[con_var]
        elif con_var in self.con_surf_to_dindex:
            avl_con_var = self.con_surf_to_dindex[con_var]
        elif con_var in self.dindex_to_con_surf:
            avl_con_var = con_var
        else:
            raise ValueError(
                f"specified contraint variable `{con_var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.con_surf_to_dindex.items()]}."
            )

        # check that the type of val is correct
        if not isinstance(val, (int, float, np.floating, np.integer)):
            raise TypeError(f"contraint value must be a int or float for contraint {var}. Got {type(val)}")

        self.avl.conset(avl_var, f"{avl_con_var} {val} \n")

    def set_trim_condition(self, variable: str, val: float):
        """Set a variable of the trim condition (analogus to the AVL's C1 command from the OPER menu)

        Args:
            variable: variable to be set. Options are ["bankAng", "CL", "velocity", "mass", "dens", "G", "X cg","Y cg","Z cg"]
            val: value to set the variable to

        """

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

        if variable not in options:
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

    def get_body_forces(self) -> Dict[str, float]:
        """Get the aerodynamic data for the last run case and return it as a dictionary.

        Returns:
            Dict[str, float]: Dictionary of aerodynamic data. The keys the aerodyanmic coefficients.
        """

        body_data = {}
        body_names = self.get_body_names()
        
        for body_name in body_names:
            body_data[body_name] = {}
        
        for key, avl_key in self.body_forces_var_to_fort_var.items():
            val = self.get_avl_fort_arr(avl_key[0], avl_key[1])
            for idx_body, body_name in enumerate(body_names):
                if (len(avl_key) > 2):
                    body_data[body_name][key] = val[idx_body,avl_key[2]]
                else:
                    body_data[body_name][key] = val[idx_body]

        return body_data

    def get_control_stab_derivs(self) -> Dict[str, float]:
        """Get the control surface derivative data, i.e. dCL/dElevator,
        for the current analysis run

        Returns:
            stab_deriv_dict: The dictionary of control surface derivatives, d{force coefficent}/d{control surface}.
        """

        deriv_data = {}

        control_names = self.get_control_names()

        for key, avl_key in self.case_derivs_to_fort_var.items():
            slicer = slice(0, len(control_names))
            val_arr = self.get_avl_fort_arr(*avl_key, slicer=slicer)
            for idx_control, val in enumerate(val_arr):
                control = control_names[idx_control]
                deriv_data[self._get_deriv_key(control, key)] = val[()]

        return deriv_data

    def get_stab_derivs(self) -> Dict[str, float]:
        """gets the stability derivates after an analysis run

        Returns:
            stab_deriv_dict: Dictionary of stability derivatives.
        """
        deriv_data = {}

        for func_key, avl_key in self.case_stab_derivs_to_fort_var.items():
            val_arr = self.get_avl_fort_arr(*avl_key)
            deriv_data[func_key] = val_arr[()]

        return deriv_data

    def get_body_axis_derivs(self) -> Dict[str, float]:
        """gets the body-axis derivates after an analysis run

        Returns:
            body_deriv_dict: Dictionary of stability derivatives.
        """
        deriv_data = {}

        for func_key, avl_key in self.case_body_derivs_to_fort_var.items():
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

    def get_avl_fort_arr(self, common_block: str, variable: str, slicer: Optional[slice] = None) -> np.ndarray:
        """Get data from the Fortran level common block data structure. see AVL.INC for all availible variables

        Args:
            common_block: Name of the common block of the variable like `CASE_R`
            variable: Name of the variable to retrive
            slicer: slice applied to the common block variable to return a subset of the data. i.e. (100) or slice(2, 5)

        Returns:
            val: value of variable after applying the slice (if present)
        """
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

        # if the array is of size 1 then just return the float
        if val.size == 1 and val.shape == ():
            val = val.flatten()[0]

        return val

    def set_avl_fort_arr(self, common_block: str, variable: str, val: float, slicer: Optional[slice] = None) -> None:
        """Set data from the Fortran level common block data structure. see AVL.INC for all availible variables

        Args:
            common_block: Name of the common block of the variable like `CASE_R`
            variable: Name of the variable to retrive
            val: value to set, which can be a numpy array
            slicer: slice applied to the common block variable to return a subset of the data. i.e. (100) or slice(2, 5)

        """
        # convert from c ordering to fortran ordering
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
            if isinstance(slicer, int):
                if isinstance(val, np.ndarray):
                    if val.size == 1:
                        # convert it to a float
                        fort_val = val.flatten()[0]
                    else:
                        raise ValueError(f"slicer {slicer} is integer, but value {val} isn't size 1")
                else:
                    fort_val = val
            elif isinstance(slicer, slice):
                # there is no need to flip the ordering for a 1D slice
                fort_val = val
            else:
                # the slicer is multidimensional
                # flip the order of the slicer to match the cordinates of the val
                slicer = slicer[::-1]
                fort_val = val

            original_val = getattr(common_block_obj, variable.upper())
            original_val[slicer] = fort_val
            setattr(common_block_obj, variable.upper(), original_val)

        return

    def get_surface_forces(self) -> Dict[str, Dict[str, float]]:
        """Returns the force data from each surface (including mirriored surfaces)

        Returns:
            surf_data_dict: a dictionary of surface data where the first key is the surface and the second is the force coefficient
        """

        # add a dictionary for each surface that will be filled later
        surf_data = {}
        for surf in self.surface_names:
            surf_data[surf] = {}

        for key, avl_key in self.case_surf_var_to_fort_var.items():
            vals = self.get_avl_fort_arr(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(self.surface_names):
                surf_data[surf_name][key] = vals[idx_surf]

        return surf_data

    def get_parameter(self, param_key: str) -> float:
        """
        Analogous to ruinont Modify parameters for the OPER menu to view parameters.

        Args:
            param_key: the name of the parameter to return

        Returns:
            param_val: the value of the parameter

        """
        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")

        # the key could be in one of two dicts
        if param_key in self.param_idx_dict:
            idx_param = self.param_idx_dict[param_key]
        elif param_key in self.state_param_idx_dict:
            idx_param = self.state_param_idx_dict[param_key]
        else:
            raise ValueError(
                f"param '{param_key}' not in possilbe list\n"
                f"{[k for k in self.param_idx_dict] + [k for k in self.state_param_idx_dict]}"
            )
        # [0] because optvl only supports 1 run case
        param_val = parvals[0][idx_param]

        return param_val

    def get_constraint(self, con_key: str) -> float:
        """Get the value of a constraint

        Args:
            con_key: name of the constraint. Options are ["alpha","beta","roll rate","pitch rate","yaw rate","CL","CY","Cl","Cm","Cl'"]

        Returns:
            con_val: value of the constraint
        """
        convals = self.get_avl_fort_arr("CASE_R", "CONVAL")

        # [0] because optvl only supports 1 run case
        con_val = convals[0][self.conval_idx_dict[con_key]]

        return con_val

    def set_parameter(self, param_key: str, param_val: float) -> None:
        """Modify a parameter of the run (analogous to M from the OPER menu in AVL).

        Args:
            param_key: parameter to modify. Options are ["CD0","bank","elevation","heading","Mach","velocity","density","grav.acc.","turn rad.","load fac.","X cg","Y cg","Z cg","mass","Ixx","Iyy","Izz","Ixy","Iyz","Izx","visc CL_a","visc CL_u","visc Cm_a","visc Cm_u"]
            param_val: value to set

        """

        # warn the user that alpha, beta,
        if param_key in ["alpha", "beta", "pb/2V", "qc/2V", "rb/2V", "CL", "roll rate", "picth rate", "yaw rate"]:
            raise ValueError(
                "alpha, beta, pb/2V, qc/2V, rb/2V, and CL are not allowed to be set,\n\
                             they are calculated during each run based on the constraints. to specify\n\
                             one of these values use the add_constraint method."
            )

        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")
        # [0] because optvl only supports 1 run case
        parvals[0][self.param_idx_dict[param_key]] = param_val

        self.set_avl_fort_arr("CASE_R", "PARVAL", parvals)

        # (1) here because we want to set the first runcase with fortran indexing (the only one)
        self.avl.set_params(1)

    def get_control_deflections(self) -> Dict[str, float]:
        """Get the deflections of all the control surfaces

        Returns:
            def_dict: dictionary of control surfaces as the keys and deflections as the values
        """
        control_surfaces = self.get_control_names()

        def_arr = copy.deepcopy(self.get_avl_fort_arr("CASE_R", "DELCON"))

        def_dict = {}
        for idx_con, con_surf in enumerate(control_surfaces):
            def_dict[con_surf] = def_arr[idx_con]

        return def_dict

    def get_control_deflection(self, con_surf) -> Dict[str, float]:
        """get the deflections of the control surfaces

        Returns:
            val: deflection of control surface
        """
        control_surfaces = self.get_control_names()
        idx_surf = control_surfaces.index(con_surf)

        def_arr = copy.deepcopy(self.get_avl_fort_arr("CASE_R", "DELCON"))

        val = def_arr[idx_surf]

        return val

    def set_control_deflection(self, con_surf, val) -> Dict[str, float]:
        """set the deflections of the control surfaces

        args:
            con_surf : Name or D value (D1, D2, etc.) of the control surface
        """
        con_surf_names = list(self.con_surf_to_dindex.keys())
        con_surf_dindex = list(self.dindex_to_con_surf.keys())

        if con_surf in con_surf_names:
            idx_con_surf = con_surf_names.index(con_surf)
        elif con_surf in con_surf_dindex:
            idx_con_surf = con_surf_dindex.index(con_surf)
            con_surf = self.dindex_to_con_surf[con_surf]

        vars_idx = self.conval_idx_dict[con_surf]

        self.set_avl_fort_arr("CASE_R", "CONVAL", val, (vars_idx,))
        self.set_avl_fort_arr("CASE_R", "DELCON", val, slicer=(idx_con_surf,))

    def set_control_deflections(self, def_dict: Dict[str, float]):
        """get the deflections of all the control surfaces

        args:
            def_dict: dictionary of control surfaces as the keys and deflections as the values
        """
        control_surfaces = self.get_control_names()

        for con in def_dict:
            idx_con = control_surfaces.index(con)
            self.set_avl_fort_arr("CASE_R", "DELCON", def_dict[con], slicer=idx_con)

    def get_hinge_moments(self) -> Dict[str, float]:
        """Get the hinge moments from the fortran layer and return them as a dictionary

        Returns:
            hinge_moments: array of control surface moments. The order the control surfaces are declared are the indices,
        """
        hinge_moments = {}

        control_surfaces = self.get_control_names()
        mom_array = self.get_avl_fort_arr("CASE_R", "CHINGE")

        for idx_con, con_surf in enumerate(control_surfaces):
            hinge_moments[con_surf] = mom_array[idx_con]

        return hinge_moments

    def get_strip_forces(self) -> Dict[str, Dict[str, np.ndarray]]:
        """Get force data for each strip (chordwise segment) of the mesh.

        Returns:
            strip_data: dictionary of strip data. The keys are ["chord", "width", "X LE", "Y LE", "Z LE", "twist","CL", "CD", "CDv", "downwash", "CX", "CY", "CZ","Cm", "Cn", "Cl","CL strip", "CD strip", "CF strip", "Cm strip","CL perp","Cm c/4,"Cm LE"]

        """
        # fmt: off
        var_to_fort_var = {
            # geometric quantities
            "X LE": ["STRP_R", "RLE", (slice(None), 0)],  # control point leading edge coordinates
            "Y LE": ["STRP_R", "RLE", (slice(None), 1)],  # control point leading edge coordinates
            "Z LE": ["STRP_R", "RLE", (slice(None), 2)],  # control point leading edge coordinates
            "chord": ["STRP_R", "CHORD"],
            "width": ["STRP_R", "WSTRIP"],
            "twist": ["STRP_R", "AINC"],

            # strip contributions to total lift and drag from strip integration
            "CL": ["STRP_R", "CLSTRP"],
            "CD": ["STRP_R", "CDSTRP"],
            "CDv" : ["STRP_R","CDV_LSTRP"],  # strip viscous drag in stability axes
            "downwash" : ["STRP_R","DWWAKE"],
            
            
            # strip contributions to non-dimensionalized forces 
            "CX": ["STRP_R", "CFSTRP", (slice(None), 0)], 
            "CY": ["STRP_R", "CFSTRP", (slice(None), 1)], 
            "CZ": ["STRP_R", "CFSTRP", (slice(None), 2)],   
            
            # strip contributions to total moments (body frame)
            "Cl": ["STRP_R", "CMSTRP", (slice(None), 0)], # previously CR
            "Cm": ["STRP_R", "CMSTRP", (slice(None), 1)], # previously CM
            "Cn": ["STRP_R", "CMSTRP", (slice(None), 2)], # previously CN
            
            
            # forces non-dimentionalized by strip quantities
            "CL strip" : ["STRP_R", "CL_LSTRP"],
            "CD strip" : ["STRP_R", "CD_LSTRP"],
            "CF strip" : ["STRP_R", "CF_LSTRP"], # forces in 3 directions
            "Cm strip" : ["STRP_R", "CM_LSTRP"], # moments in 3 directions

            # additional forces and moments
            "CL perp" : ["STRP_R", "CLT_LSTRP"], # strip CL referenced to Vperp,
            "Cm c/4" : ["STRP_R","CMC4_LSTRP"],  # strip pitching moment about c/4 and
            "Cm LE" : ["STRP_R","CMLE_LSTRP"],  # strip pitching moment about LE vector
            "spanloading" : ["STRP_R","CNC"],   # strip spanloading 
            
            # TODO: add
            #  & CF_LSTRP(3,NSMAX),   CM_LSTRP(3,NSMAX),    ! strip forces in body axes referenced to strip area and 1/4 chord
        }    
        # fmt: on

        # add a dictionary for each surface that will be filled later
        strip_data = {}
        for surf in self.surface_names:
            strip_data[surf] = {}

        for key, avl_key in var_to_fort_var.items():
            vals = self.get_avl_fort_arr(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(self.surface_names):
                idx_srp_beg, idx_srp_end = self._get_surface_strip_indices(idx_surf)
                strip_data[surf_name][key] = vals[idx_srp_beg:idx_srp_end]

        # process the data
        # add a few more pieces of data that are output on the fortran layer
        for surf_key in strip_data:
            # convert the twist to degrees
            strip_data[surf_key]["twist"] = 180 / np.pi * strip_data[surf_key]["twist"]
            # add the area of each strip
            strip_data[surf_key]["area"] = strip_data[surf_key]["width"] * strip_data[surf_key]["chord"]

            # formula is directly from AVL
            xcp = 0.25 - strip_data[surf_key]["Cm c/4"] / strip_data[surf_key]["CL perp"]
            strip_data[surf_key]["CP x/c"] = xcp

        # get length of along the surface of each strip
        for idx_surf, surf_key in enumerate(strip_data):
            xles = strip_data[surf_key]["X LE"]
            yles = strip_data[surf_key]["Y LE"]
            zles = strip_data[surf_key]["Z LE"]

            n_strips = len(yles)
            sles = np.zeros(n_strips)

            sles[0] = 0
            for idx_strip in range(1, n_strips):
                dx = xles[idx_strip] - xles[idx_strip - 1]
                dy = yles[idx_strip] - yles[idx_strip - 1]
                dz = zles[idx_strip] - zles[idx_strip - 1]

                sles[idx_strip] = sles[idx_strip - 1] + np.sqrt(dx**2 + dy**2 + dz**2)

            strip_data[surf_key]["S LE"] = sles

        ref_data = self.get_reference_data()
        cref = ref_data["Cref"]
        bref = ref_data["Bref"]

        for surf_key in strip_data:
            # add sectional lift and drag
            strip_data[surf_key]["lift dist"] = strip_data[surf_key]["CL"] * strip_data[surf_key]["chord"] / cref
            strip_data[surf_key]["drag dist"] = strip_data[surf_key]["CD"] * strip_data[surf_key]["chord"] / cref
            strip_data[surf_key]["roll dist"] = strip_data[surf_key]["Cl"] * (strip_data[surf_key]["chord"] / cref) ** 2
            strip_data[surf_key]["yaw dist"] = (
                strip_data[surf_key]["Cn"] * strip_data[surf_key]["chord"] ** 2 / (bref * cref)
            )

        return strip_data

    def _get_surface_strip_indices(self, idx_surf: int):
        num_strips = np.trim_zeros(self.get_avl_fort_arr("SURF_I", "NJ"))
        idx_srp_beg = np.sum(num_strips[:idx_surf])
        idx_srp_end = np.sum(num_strips[: idx_surf + 1])

        return idx_srp_beg, idx_srp_end

    # region --- modal analysis api
    def execute_eigen_mode_calc(self):
        """Execute a modal analysis (x from the MODE menu in AVL)"""
        self.avl.execute_eigenmode_calc()

    def get_eigenvalues(self) -> np.ndarray:
        """After running an eigenmode calculation, this function will return the eigenvalues in the order used by AVL

        Returns:
            eig_vals: array of eigen values
        """

        # get the number of "valid" eigenvalues from avl
        # [0] because optvl only supports 1 run case
        num_eigen = self.get_avl_fort_arr("CASE_I", "NEIGEN")[0]

        # 0 because optvl only supports 1 run case
        slicer = (0, slice(0, num_eigen))
        # get the eigenvalues from avl
        eig_vals = self.get_avl_fort_arr("CASE_Z", "EVAL", slicer=slicer)
        return eig_vals

    def get_eigenvectors(self) -> np.ndarray:
        """After running an eigenmode calculation, this function will return the eigenvalues in the order used by AVL

        Returns:
            eig_vec: 2D array of eigen vectors
        """

        # get the number of "valid" eigenvalues from avl
        # [0] because optvl only supports 1 run case
        num_eigen = self.get_avl_fort_arr("CASE_I", "NEIGEN")[0]

        # 0 because optvl only supports 1 run case
        slicer = (0, slice(0, num_eigen), slice(None))
        eig_vecs = self.get_avl_fort_arr("CASE_Z", "EVEC", slicer=slicer)

        return eig_vecs

    def get_system_matrix(self, in_body_axis=False) -> np.ndarray:
        """ returns just the A system matrix used for the eigenmode calculation"""
        # this routine is mostly to not introduce breaking chnages
        
        # call the routine for getting them all, but just ignore the B and R parts
        Asys, _, _ = self.get_system_matrices(in_body_axis=in_body_axis)
        
        return Asys

    def get_system_matrices(self, in_body_axis=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """returns the A,B, and R system matrices used in amode.f
        
        args:
            in_body_axis: apply the sign changes to the matrix to put it in the body axis
        
        Returns:
            Asys: 2D array representing the system matrix for the eigen value analysis
            Bsys: 2D array representing the system matrix for control surfaces
            Rsys: 1D array representing the RHS of the dynamics equation
        """
        
        # get the dimesion of the A matrix from the eig_vals    
        eig_vals = self.get_avl_fort_arr("CASE_Z", "EVAL")
        jemax = eig_vals.shape[1]
        Asys = np.zeros((jemax,jemax), order="F")
        Bsys = np.zeros((jemax,self.NDMAX), order="F")
        Rsys = np.zeros((jemax), order="F")
        
        # 1 because optvl only supports 1 run case and we are using fortran base 1 indexing
        irun_case = 1
        self.avl.get_system_matrices(irun_case,Asys, Bsys, Rsys)
        
        num_controls = self.get_num_control_surfs()
        # trim the columns of the Bsys matrix
        Bsys = Bsys[:, 0:num_controls]
        
        def apply_state_signs(Asys, Bsys, Rsys):
            """
            Apply sign changes to the state matrix A and return the modified version.
            """
            nsys = Asys.shape[0]

            # Indices for sign flip
            jeu   = 0
            jew   = 1
            jeq   = 2
            jeth  = 3
            jev   = 4
            jep   = 5
            jer   = 6
            jeph  = 7
            jex   = 8
            jey   = 9
            jez   = 10
            jeps  = 11

            usgn = np.ones(nsys)
            for idx in [jeu, jew, jep, jer, jex, jez]:
                if 0 <= idx < nsys:
                    usgn[idx] = -1.0

            # Allocate result
            Asys_signed = np.zeros_like(Asys)
            Bsys_signed = np.zeros_like(Bsys)
            Rsys_signed = np.zeros_like(Rsys)

            # Apply row/column sign flips explicitly
            for i in range(nsys):
                for j in range(nsys):
                    Asys_signed[i, j] = Asys[i, j] * usgn[i] * usgn[j]
                
                for j in range(num_controls):
                    Bsys_signed[i,j] = Bsys[i, j] * usgn[i]
                
                Rsys_signed = Rsys[i]*usgn[i]

            return Asys_signed, Bsys_signed, Rsys_signed
        
        if in_body_axis:
            Asys, Bsys, Rsys = apply_state_signs(Asys, Bsys, Rsys)
        
        return Asys, Bsys, Rsys

    # region --- geometry api
    def get_control_names(self) -> List[str]:
        """Get the names of the control surfaces

        Returns:
            control_names: list of control surface names
        """
        fort_names = self.get_avl_fort_arr("CASE_C", "DNAME")
        control_names = self._fort_char_array_to_str_list(fort_names)
        return control_names

    def get_design_var_names(self) -> List[str]:
        """Get the names of the design_var surfaces

        Returns:
            design_var_names: list of design_var surface names
        """
        fort_names = self.get_avl_fort_arr("CASE_C", "GNAME")
        design_var_names = self._fort_char_array_to_str_list(fort_names)
        return design_var_names

    def get_surface_names(self, remove_dublicated: Optional[bool] = False) -> List[str]:
        """Get the surface names from the geometry

        Args:
            remove_dublicated: remove the surface that were created by duplication about symmetry planes

        Returns:
            surf_names: list of surface names
        """
        fort_names = self.get_avl_fort_arr("CASE_C", "STITLE")
        surf_names = self._fort_char_array_to_str_list(fort_names)

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

    def get_body_names(self, remove_dublicated: Optional[bool] = False) -> List[str]:
        """Get the body names from the geometry

        Args:
            remove_dublicated: remove the body that were created by duplication about symmetry planes

        Returns:
            body_names: list of body names
        """
        fort_names = self.get_avl_fort_arr("CASE_C", "BTITLE")
        body_names = self._fort_char_array_to_str_list(fort_names)

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

    def get_con_surf_param(self, surf_name: str, idx_sec: int, param: str) -> np.ndarray:
        """Returns the parameters that define the control surface. Can also get design variables (AVL).

        Args:
            surf_name: the name of the surface containing the control surface
            idx_sec: the section index of the control surface data
            param: control surface parameter to get

        Returns:
            parm: parameter value
        """
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.con_surf_to_fort_var[surf_name].keys():
            fort_var = self.con_surf_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has control surface data {self.con_surf_to_fort_var[surf_name].keys()}"
            )
        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_sec])
        return param

    def __get_des_var_param(self, surf_name: str, idx_sec: int, param: str) -> np.ndarray:
        """Returns the parameters that define the control surface. Can also get design variables (AVL).

        Args:
            surf_name: the name of the surface containing the control surface
            idx_sec: the section index of the control surface data
            param: control surface parameter to get

        Returns:
            parm: parameter value
        """
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.des_var_to_fort_var[surf_name].keys():
            fort_var = self.des_var_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has design variable data {self.des_var_to_fort_var[surf_name].keys()}"
            )
        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_sec])
        return param

    def set_con_surf_param(
        self, surf_name: str, idx_slice: int, param: str, val: float, update_geom: Optional[bool] = True
    ):
        """Sets the parameters that define the control surface. Can also set design variables (AVL).

        Args:
            surf_name: the name of the surface containing the control surface
            idx_slice: the section index of the control surface data
            param: control surface parameter to set
            val: value to set
            update_geom: flag to update the geometry after setting

        """
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.con_surf_to_fort_var[surf_name].keys():
            fort_var = self.con_surf_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has control surface data {self.con_surf_to_fort_var[surf_name].keys()}"
            )
        # param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_slice])
        self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2][idx_slice])

        if update_geom:
            self.avl.update_surfaces()

    def get_surface_param(self, surf_name: str, param: str) -> np.ndarray:
        """Get a parameter of a specified surface. Does not get control surface or design variables.

        Args:
            surf_name: the surface containing the parameter
            param: the surface parameter to return. Could be either geometric or paneling


        Returns:
            param: the parameter of the surface
        """

        # check that param is in self.surf_geom_to_fort_var
        if param in self.surf_section_geom_to_fort_var[surf_name].keys():
            # TODO-JLA: handle section parameters differently
            # return the data for each section
            fort_vars = self.surf_section_geom_to_fort_var[surf_name][param]
            param_list = []

            for slicer in fort_vars[2]:
                param = self.get_avl_fort_arr(fort_vars[0], fort_vars[1], slicer=slicer)

                param_list.append(copy.deepcopy(param))

            return param_list

        if param in self.surf_geom_to_fort_var[surf_name].keys():
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
        elif param in self.surf_pannel_to_fort_var[surf_name].keys():
            fort_var = self.surf_pannel_to_fort_var[surf_name][param]
        elif param in self.con_surf_to_fort_var[surf_name].keys():
            warnings.warn(
                "OptVL WARNING - Getting control surface and design variables is not supported with this function.\n"
                "Use the get_con_surf_param function.",
                stacklevel=2,
            )
        elif param in ["afiles", "airfoils", "naca", "xfminmax"]:
            warnings.warn(
                "OptVL WARNING - Getting section geometry data using airfoils, afiles, or naca with this function is not supported.",
                stacklevel=2,
            )
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_section_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )

        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2])
        return copy.deepcopy(param)  # return the value of the array, but not a reference to avoid sideffects

    def set_surface_param(self, surf_name: str, param: str, val: float, update_geom: bool = True):
        """Set a parameter of a specified surface. Supports setting params related to geometry and panelling.
        Section geometry can be directly set here but this is not recommended. Use set_section_coordinates instead.

        Args:
            surf_name: the surface containing the parameter
            param: the surface parameter to return. Could be either geometric or paneling
            val: value to set
            update_geom: flag to update the geometry after setting
        """

        # check that the surface is in the set of unique surfaces and not duplicated
        if surf_name not in self.unique_surface_names:
            raise ValueError(
                f"Only non-duplicates surface parameters can be set, {surf_name} not found in {self.unique_surface_names}"
            )

        # check that param is in self.surf_geom_to_fort_var
        if param in self.surf_section_geom_to_fort_var[surf_name].keys():
            # TODO-JLA: handle section parameters differently
            # return the data for each section
            # Set surface cross section geometry variables (not recommended)
            warnings.warn(
                "OptVL WARNING - Updating section geometry data directly is not recommened!\n "
                "OptVL will not verify that the camber line is consistent with the given coordiantes.\n"
                "Use the set_section_coordinates function.",
                stacklevel=2,
            )

            fort_vars = self.surf_section_geom_to_fort_var[surf_name][param]

            for idx_sec, slicer in enumerate(fort_vars[2]):
                self.set_avl_fort_arr(fort_vars[0], fort_vars[1], val[idx_sec], slicer=slicer)

        elif param in self.surf_geom_to_fort_var[surf_name].keys():
            # Set basic surface geometry variables
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
            self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2])

        elif param in self.surf_pannel_to_fort_var[surf_name].keys():
            # Set surface panelling variables
            fort_var = self.surf_pannel_to_fort_var[surf_name][param]
            self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2])
        elif param in ["afiles", "airfoils", "naca", "xfminmax"]:
            # Cannot indirectly update the cross sections like this. Would over complicate this function when it can easily be handled by set_section_coordinates
            warnings.warn(
                "OptVL WARNING - Updating section geometry data with using airfoils, afiles, or naca is not supported with this function.\n"
                "Use the set_section_coordinates function.",
                stacklevel=2,
            )
            pass
        elif param in self.con_surf_to_fort_var[surf_name].keys():
            # the control surface and design variables need to be handeled differently because the number at each section is variable
            warnings.warn(
                "OptVL WARNING - Updating control surface and design variables is not supported with this function.\n"
                "Use the set_con_surf_param function.",
                stacklevel=2,
            )
            pass
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )

        if update_geom:
            self.avl.update_surfaces()

    def get_surface_params(
        self,
        include_geom: bool = True,
        include_section_geom: bool = False,
        include_paneling: bool = False,
        include_con_surf: bool = False,
        include_des_vars: bool = False,
        include_airfoils: bool = False,
    ) -> Dict[str, Dict[str, Any]]:
        """Get all the surface level parameters for each suface

        Args:
            include_geom: flag to include geometry data in the output. The data is ["scale", "translate", "angle", "xles", "yles", "zles", "chords", "aincs", "clcdsec", "claf"]
            include_section_geom: flag to include section geometry data in the output. The data is ["xasec", "sasec", "tasec", xuasec, xlasec, zlasec, zuasec, casec, nasec]
            include_paneling:  flag to include paneling information in the output. The data ["nchordwise", "cspace","nspan", "sspace","sspaces","nspans","yduplicate", "wake", "able", "load", "use surface spacing", "component"]
            include_con_surf:  flag to include control surface and design variable data in the output. This is data like the hinge vector and gain.
            include_airfoils:  flag to include airfoil file data in the output

        Return:
            surf_data: Nested dictionary where the 1st key is the surface name and the 2nd key is the parameter.
        """
        surf_data = {}

        for surf_name in self.unique_surface_names:
            surf_data[surf_name] = {}
            if include_geom:
                for var in self.surf_geom_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

            if include_section_geom:
                # add section geometry parameters if requested
                for var in self.surf_section_geom_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

            idx_surf = self.surface_names.index(surf_name)
            if include_paneling:
                # add paneling parameters if requested
                for var in self.surf_pannel_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

                if not self.get_avl_fort_arr("SURF_GEOM_L", "LDUPL")[idx_surf]:
                    surf_data[surf_name].pop("yduplicate")

            if include_con_surf:
                # add control surface and design variable parameters if requested
                for var in self.con_surf_to_fort_var[surf_name]:
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]

                    slice_data = []
                    for idx_sec in range(num_sec):
                        tmp = self.get_con_surf_param(surf_name, idx_sec, var)
                        slice_data.append(tmp)

                    surf_data[surf_name][var] = slice_data

            if include_des_vars:
                # add control surface and design variable parameters if requested
                for var in self.des_var_to_fort_var[surf_name]:
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]

                    slice_data = []
                    for idx_sec in range(num_sec):
                        tmp = self.__get_des_var_param(surf_name, idx_sec, var)
                        slice_data.append(tmp)

                    surf_data[surf_name][var] = slice_data

            if include_airfoils:
                # add airfoil files/coordinates/naca names if requested
                # NOTE: If section geometry was modified NACA and afile data won't be accurate
                # NOTE: If NACA is set then airfoil will be zeros. Check the raw section data instead
                afiles = []
                nacas = []
                airfoils = []
                xfminmax = []

                num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]
                for idx_sec in range(num_sec):
                    afile = self.__fort_char_array_to_str(self.avl.CASE_C.AFILES[idx_sec, idx_surf])
                    afiles.append(afile)

                    naca = self.__fort_char_array_to_str(self.avl.CASE_C.NACA[idx_sec, idx_surf])
                    nacas.append(naca)

                    airfoilx = np.trim_zeros(
                        self.get_avl_fort_arr("SURF_GEOM_R", "XSEC")[idx_surf, idx_sec, slice(None)]
                    )
                    airfoily = np.trim_zeros(
                        self.get_avl_fort_arr("SURF_GEOM_R", "YSEC")[idx_surf, idx_sec, slice(None)]
                    )

                    airfoils.append(np.array([airfoilx, airfoily]))

                    xfmin = self.get_avl_fort_arr("SURF_GEOM_R", "XFMIN_R")[idx_surf, idx_sec]
                    xfmax = self.get_avl_fort_arr("SURF_GEOM_R", "XFMAX_R")[idx_surf, idx_sec]

                    xfminmax.append([xfmin, xfmax])

                surf_data[surf_name]["naca"] = nacas
                
                surf_data[surf_name]["airfoils"] = airfoils
                surf_data[surf_name]["afiles"] = afiles
                surf_data[surf_name]["xfminmax"] = np.array(xfminmax)

        return surf_data

    def set_surface_params(self, surf_data: Dict[str, Dict[str, any]]):
        """Set the given data of the current geometry.
        ASSUMES THE CONTROL SURFACE DATA STAYS AT THE SAME LOCATION
        (i.e  you didn't move the control surfaces to new sections or surfaces. If so re-initialize OptVL)

        Args:
            surf_data: Nested dictionary where the 1st key is the surface name and the 2nd key is the parameter.

        """
        for surf_name in surf_data:
            if surf_name not in self.unique_surface_names:
                raise ValueError(
                    f"""surface name, {surf_name}, not found in the current avl object."
                        Note duplicated surfaces can not be set directly.
                        Surfaces in file {self.unique_surface_names}
                        {self.surface_names}"""
                )

            for var in surf_data[surf_name]:
                # do not set the data this way if it is a control surface or airfoil
                if var not in [self.con_surf_to_fort_var[surf_name], "airfoils", "afiles", "naca"]:
                    self.set_surface_param(surf_name, var, surf_data[surf_name][var], update_geom=False)
                elif var in [self.con_surf_to_fort_var[surf_name]]:
                    idx_surf = self.surface_names.index(surf_name)
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]
                    slice_data = []
                    for idx_sec in range(num_sec):
                        self.set_con_surf_param(
                            surf_name, idx_sec, var, surf_data[surf_name][var][idx_sec], update_geom=False
                        )
                elif var in ["airfoils", "afiles", "naca"]:
                    # Cannot indirectly update the cross sections like this. Would over complicate this function when it can easily be handled by set_section_coordinates
                    warnings.warn(
                        "OptVL WARNING - Updating section geometry data with using airfoils, afiles, or naca is not supported with this function.\n"
                        "Use the set_section_coordinates function.",
                        stacklevel=2,
                    )
                else:
                    pass

        # update the geometry once at the end
        self.avl.update_surfaces()

    def get_body_param(self, body_name: str, param: str) -> np.ndarray:
        """Get a parameter of a specified body

        Args:
            body_name: the body containing the parameter
            param: the body parameter to return. Could be either geometric or paneling


        Returns:
            val: the val of parameter of the body
        """

        body_names = self.get_body_names()
        idx_body = body_names.index(body_name)
        unique_body_names = self.get_body_names(remove_dublicated=True)
        if body_name not in unique_body_names:
            raise ValueError(
                f"""body name, {body_name}, not found in the current avl object."
                    Note duplicated bodies can not be set directly.
                    Surfaces in file {unique_body_names}
                    {body_names}"""
            )

        # yduplicate is the only parameter that may not exist for a body, everything else is required
        if not self.get_avl_fort_arr(
            self.body_geom_to_fort_var["yduplicate"][0], self.body_geom_to_fort_var["yduplicate"][1]
        )[idx_body]:
            raise ValueError(f"{param}, not set for {body_name}")

        val = self.get_avl_fort_arr(
            self.body_geom_to_fort_var[param][0], self.body_geom_to_fort_var[param][1], slicer=idx_body
        )

        if param == "bfile":
            val = self.__fort_char_array_to_str(val)

        return copy.deepcopy(val)  # return the value of the array, but not a reference to avoid sideffects

    def set_body_param(self, body_name: str, param: str, val, update_geom: bool = True):
        """Set a parameter of a specified body

        Args:
            body_name: the body containing the parameter
            param: the surface parameter to return. Could be either geometric or paneling
            val: value to set
            update_geom: flag to update the geometry after setting
        """

        body_names = self.get_body_names()
        idx_body = body_names.index(body_name)
        unique_body_names = self.get_body_names(remove_dublicated=True)
        if body_name not in unique_body_names:
            raise ValueError(
                f"""body name, {body_name}, not found in the current avl object."
                    Note duplicated bodies can not be set directly.
                    Surfaces in file {unique_body_names}
                    {body_names}"""
            )

        if param in self.body_geom_to_fort_var.keys():
            fort_var = self.body_geom_to_fort_var[param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {body_name}, that has geom data {list(self.body_geom_to_fort_var.keys())}"
            )

        # Set the slicer
        if isinstance(val, np.ndarray):
            slicer = (idx_body, slice(0, len(val)))
        else:
            slicer = idx_body

        # Handle strings differently
        if isinstance(val, str):
            operator.setitem(getattr(getattr(self.avl, fort_var[0]), fort_var[1]), idx_body, val)
        else:
            self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=slicer)

        if update_geom:
            self.avl.update_bodies()

    def get_body_params(
        self,
        include_body_oml: bool = False,
    ) -> Dict[str, Dict[str, Any]]:
        """Get the parameters of the bodies

        Args:
            include_body_oml: include the raw oml coordinates in the output dict

        Returns:
            body_data: Nested dictionary where the 1st key is the body name and the 2nd key is the parameter.
        """

        body_names = self.get_body_names()
        unique_body_names = self.get_body_names(remove_dublicated=True)
        body_data = {}

        for body_name in unique_body_names:
            idx_body = body_names.index(body_name)
            body_data[body_name] = {}

            for var, fort_var in self.body_geom_to_fort_var.items():
                val = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=idx_body)  #
                if var == "bfile":
                    val = self.__fort_char_array_to_str(val)

                body_data[body_name][var] = val

            if not self.get_avl_fort_arr(
                self.body_geom_to_fort_var["yduplicate"][0], self.body_geom_to_fort_var["yduplicate"][1]
            )[idx_body]:
                body_data[body_name].pop("yduplicate")

            if include_body_oml:
                body_data[body_name]["body_oml"] = np.array(
                    [
                        np.trim_zeros(self.get_avl_fort_arr("BODY_GEOM_R", "XBOD_R")[idx_body, slice(None)]),
                        np.trim_zeros(self.get_avl_fort_arr("BODY_GEOM_R", "YBOD_R")[idx_body, slice(None)]),
                    ]
                )

        return body_data

    def set_body_params(self, body_data: Dict[str, Dict[str, Any]]):
        """Set the give body data of the current geometry.

        Args:
            body_data: Nested dictionary where the 1st key is the body name and the 2nd key is the parameter.

        """
        body_names = self.get_body_names()
        unique_body_names = self.get_body_names(remove_dublicated=True)
        for body_name in body_data:
            if body_name not in unique_body_names:
                raise ValueError(
                    f"""body name, {body_name}, not found in the current avl object."
                        Note duplicated bodies can not be set directly.
                        Surfaces in file {unique_body_names}
                        {body_names}"""
                )

            for var in body_data[body_name]:
                if var == "body_oml":
                    raise RuntimeError("The body oml cannot be set with this function! Use set_body_coordinates.")
                self.set_body_param(body_name, var, body_data[body_name][var], update_geom=False)

        # update the geometry once at the end
        self.avl.update_bodies()

    def get_header_params(self) -> Dict:
        """Gets the header input settings in AVL and returns them in a dictionary.

        Returns:
            Dict[str]: Dictionary containing the header input settings in AVL
        """
        header_data = {}

        header_data["title"] = self.avl.CASE_C.TITLE

        for var, fort_var in self.header_var_to_fort_var.items():
            val = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=None)
            header_data[var] = val

        return header_data

    def get_input_dict(
        self, include_surfaces: bool = True, include_section_geom: bool = False, include_bodies: bool = True
    ) -> Dict[str, Dict[str, Any]]:
        """Returns all input information from AVL in input dictionary format.

        Args:
            include_surfaces (bool, optional): Include all surfaces in the dictionary. Defaults to True.
            include_section_geom (bool, optional): Include all the section geometry information for each surface. Defaults to False.
            include_bodies (bool, optional): Include all bodies in the dictionary. Defaults to True.

        Returns:
            Dict[str, Dict[str, Any]]: OptVL input dictionary
        """

        input_dict = self.get_header_params()

        if include_surfaces:
            input_dict["surfaces"] = self.get_surface_params(
                include_paneling=True,
                include_airfoils=True,
                include_con_surf=True,
                include_section_geom=include_section_geom,
            )
        else:
            input_dict["surfaces"] = {}

        if include_bodies:
            input_dict["bodies"] = self.get_body_params(include_body_oml=True)
        else:
            input_dict["bodies"] = {}

        return input_dict

    # region --- geometry file writing api
    def write_geom_file(self, filename: str):
        """Write the current geometry to a file

        Args:
            filename: name of the output AVL-style geometry file
        """
        with open(filename, "w") as fid:
            # write the header
            fid.write(f"# generated using OptVL v{__version__}\n")
            self._write_header(fid)

            surf_data = self.get_surface_params(
                include_geom=True,
                include_section_geom=True,
                include_paneling=True,
                include_con_surf=True,
                include_des_vars=True,
                include_airfoils=True,
            )

            for surf_name in surf_data:
                self._write_surface(fid, surf_name, surf_data[surf_name])

            body_data = self.get_body_params()
            for body_name in body_data:
                self._write_body(fid, body_name, body_data[body_name])

    def _write_fort_vars(self, fid: TextIO, common_block: str, fort_var: str, newline: bool = True) -> None:
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

    def _write_header(self, fid: TextIO):
        """Write the header to a file"""
        # write the name of the aircraft

        self._write_banner(fid, "Header")
        title_array = self.get_avl_fort_arr("case_c", "title")
        title = self.__fort_char_array_to_str(title_array)
        fid.write(f"{title}\n")

        fid.write("#Mach\n")
        self._write_fort_vars(fid, "case_r", "mach0")

        fid.write("#IYsym   IZsym   Zsym\n")
        self._write_fort_vars(fid, "case_i", "iysym", newline=False)
        self._write_fort_vars(fid, "case_i", "izsym", newline=False)
        self._write_fort_vars(fid, "case_r", "zsym")

        fid.write("#Sref    Cref    Bref\n")
        self._write_fort_vars(fid, "case_r", "sref", newline=False)
        self._write_fort_vars(fid, "case_r", "cref", newline=False)
        self._write_fort_vars(fid, "case_r", "bref")

        fid.write("#Xref    Yref    Zref\n")
        self._write_fort_vars(fid, "case_r", "XYZREF0")  # The inpute routine reads the defaults

        fid.write("#CD0\n")
        self._write_fort_vars(fid, "case_r", "CDREF0")

        # fid.write(f" {self.get_avl_fort_arr('case_r', 'sref')}")

    def _write_banner(self, fid: TextIO, header: str, line_width: int = 80):
        header = " " + header + " "  # pad with spaces

        width = line_width - 1

        banner = f"#{'=' * (width)}\n#{header.center(width, '-')}\n#{'=' * (width)}\n"
        fid.write(banner)

    def _write_body(self, fid: TextIO, body_name: str, data: Dict[str, float]):
        self._write_banner(fid, body_name)
        fid.write("BODY\n")
        fid.write(f"{body_name}\n")
        fid.write("#N  Bspace\n")
        fid.write(f"{data['nvb']} {data['bspace']}\n")
        if "yduplicate" in data.keys():
            fid.write("YDUPLICATE\n")
            fid.write(f"{data['yduplicate']}\n")
        fid.write("SCALE\n")
        fid.write(f"{data['scale'][0]} {data['scale'][1]} {data['scale'][2]}\n")
        fid.write("TRANSLATE\n")
        fid.write(f"{data['translate'][0]} {data['translate'][1]} {data['translate'][2]}\n")
        if data["bfile"] != "":
            fid.write("BFILE\n")
            fid.write(f"{data['bfile']}\n")

    def _write_surface(self, fid: TextIO, surf_name: str, data: Dict[str, float]):
        """Write a surface to a file"""
        # TODO add NACA and CLAF keyword support

        def _write_data(key_list: List[str], newline: bool = True):
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
        self._write_banner(fid, surf_name)
        fid.write("SURFACE\n")
        fid.write(f"{surf_name}\n")

        fid.write("#Nchordwise  Cspace  [Nspanwise  Sspace]\n")
        _write_data(["nchordwise", "cspace"], newline=False)
        if data["use surface spacing"]:
            _write_data(["nspan", "sspace"])
        else:
            fid.write("\n")

        if "yduplicate" in data:
            fid.write("YDUPLICATE\n")
            _write_data(["yduplicate"])

        if not data["wake"]:
            fid.write("NOWAKE\n")
        if not data["albe"]:
            fid.write("NOALBE\n")
        if not data["load"]:
            fid.write("NOLOAD\n")

        idx_surf = self.get_surface_index(surf_name)
        if idx_surf + 1 != data["component"]:
            # only add component keys if we have to to avoid freaking
            # people out who don't expect to see them
            fid.write("COMPONENT\n")
            _write_data(["component"])

        fid.write("SCALE\n")
        _write_data(["scale"])

        fid.write("TRANSLATE\n")
        _write_data(["translate"])

        fid.write("ANGLE\n")
        _write_data(["angle"])

        if (data["clcd"] != 0.0).any():
            fid.write(" CDCL\n")
            fid.write(
                f" {data['clcd'][0]:.6f} "
                f" {data['clcd'][1]:.6f} "
                f" {data['clcd'][2]:.6f} "
                f" {data['clcd'][3]:.6f} "
                f" {data['clcd'][4]:.6f} "
                f" {data['clcd'][5]:.6f}\n"
            )

        fid.write("#---------------------------------------\n")

        num_sec = data["chords"].size
        control_names = self.get_control_names()
        design_var_names = self.get_design_var_names()

        for idx_sec in range(num_sec):
            fid.write("SECTION\n")
            fid.write("#Xle      Yle      Zle      | Chord    Ainc     Nspan  Sspace\n")
            fid.write(
                f" {data['xles'][idx_sec]:.6f} "
                f"{data['yles'][idx_sec]:.6f} "
                f"{data['zles'][idx_sec]:.6f}   "
                f"{data['chords'][idx_sec]:.6f} "
                f"{data['aincs'][idx_sec]:.6f} "
            )

            if data["nspans"][idx_sec] != 0:
                # if no section spacing is read in by avl nspans will be zero
                fid.write(f"{data['nspans'][idx_sec]}      {data['sspaces'][idx_sec]}\n")
            elif data["use surface spacing"]:
                fid.write("\n")
            else:
                raise RuntimeError(
                    f"neither surface nor section spacing information written"
                    f"for surface {surf_name} section {idx_sec + 1}"
                )

            afile = data["afiles"][idx_sec]
            naca = data["naca"][idx_sec]
            airfoil = data["airfoils"][idx_sec]

            if naca:
                fid.write("#NACA | X1 X2\n")
                fid.write(f" NACA   {data['xfminmax'][idx_sec][0]}  {data['xfminmax'][idx_sec][1]}\n")
                fid.write(f" {naca}\n")

            if np.any(airfoil) and afile == '':
                fid.write("#AIRFOIL | X1 X2\n")
                fid.write(f"AIRFOIL   {data['xfminmax'][idx_sec][0]}  {data['xfminmax'][idx_sec][1]}\n")
                for i in range(0, min(self.IBX, len(airfoil[0, :]))):
                    fid.write(f" {airfoil[0, i]} {airfoil[1, i]}\n")

            if afile:
                fid.write("#AFILE | X1 X2\n")
                fid.write(f" AFILE {data['xfminmax'][idx_sec][0]}  {data['xfminmax'][idx_sec][1]} \n")
                fid.write(f" {afile}\n")

            # output claf and  clcd if not default
            if data["claf"][idx_sec] != 1.0:
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
                fid.write(f" {control_names[idx_cont_surf - 1]} ")
                fid.write(f" {data['gaind'][idx_sec][idx_local_cont_surf]}")
                fid.write(f" {data['xhinged'][idx_sec][idx_local_cont_surf]}")
                vhinge = data["vhinged"][idx_sec][idx_local_cont_surf]
                fid.write(f" {vhinge[0]:.6f} {vhinge[1]:.6f} {vhinge[2]:.6f}")
                fid.write(f" {data['refld'][idx_sec][idx_local_cont_surf]}\n")

            # check for design variables
            for idx_local_des_var, idx_des_var in enumerate(data["idestd"][idx_sec]):
                fid.write(" DESIGN\n")
                fid.write("#surface   gain\n")
                fid.write(f" {design_var_names[idx_des_var - 1]} ")
                fid.write(f" {data['gaing'][idx_sec][idx_local_des_var]}\n")

    # region --- Utility functions
    def get_num_surfaces(self) -> int:
        """Returns the number of surface including duplicated

        Returns:
            val: number of surfaces
        """
        """Get the number of surfaces in the geometry"""
        return self.get_avl_fort_arr("CASE_I", "NSURF")

    def get_surface_index(self, surf_name: str) -> int:
        """Given a surface name returns the index

        Args:
            surf_name: name of the surface

        Returns:
            idx_surf: index of the surface
        """
        surf_names = self.surface_names
        idx_surf = surf_names.index(surf_name)
        return idx_surf

    def get_body_index(self, body_name: str) -> int:
        """Given a body name returns the index

        Args:
            body_name: name of the body

        Returns:
            idx_body: index of the body
        """
        body_names = self.body_names
        idx_body = body_names.index(body_name)
        return idx_body

    def get_num_sections(self, surf_name: str) -> int:
        """Get the number of sections in a surface.

        Args:
            surf_name: name of the surface

        Returns:
            nsec: numer of sections
        """
        idx_surf = self.get_surface_index(surf_name)
        slice_idx_surf = (idx_surf,)
        return self.get_avl_fort_arr("SURF_GEOM_I", "NSEC", slicer=slice_idx_surf)

    def get_num_strips(self) -> int:
        """
        Get the number of strips in the mesh

        """
        return self.get_avl_fort_arr("CASE_I", "NSTRIP")

    def get_num_control_surfs(self) -> int:
        """Get the number of control surfaces

        Returns:
            val: number of control surfaces
        """
        return self.get_avl_fort_arr("CASE_I", "NCONTROL")

    def get_mesh_size(self) -> int:
        """Get the number of vortices in the mesh

        Returns:
            val: the number of vortices
        """
        return int(self.get_avl_fort_arr("CASE_I", "NVOR"))

    def get_mesh_data(self) -> Dict[str, int]:
        """Get the number of vortices in the mesh

        Returns:
            mesh_data : dict
        """

        mesh_data = {
            "bodies": self.get_avl_fort_arr("CASE_I", "NBODY"),
            "surfaces": self.get_avl_fort_arr("CASE_I", "NSURF"),
            "strips": self.get_avl_fort_arr("CASE_I", "NSTRIP"),
            "vortices": self.get_avl_fort_arr("CASE_I", "NVOR"),
            "control variables": self.get_avl_fort_arr("CASE_I", "NCONTROL"),
            "design variables": self.get_avl_fort_arr("CASE_I", "NDESIGN"),
        }

        return mesh_data

    def _str_to_fort_str(self, py_string, num_max_char):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takes a strings and returns the char array.
        """

        arr = np.zeros((), dtype=f"|S{num_max_char}")
        # pad the array with spaces
        arr = py_string + " " * (num_max_char - (len(py_string)))

        # for j in range(len(py_string)):
        #     # copy over each character
        #     arr[j] = py_string[j]

        return arr

    def _str_to_fort_char_array(self, py_string, num_max_char):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takes a strings and returns the char array.
        """

        arr = np.zeros(1, dtype=f"|S{num_max_char}")
        # pad the array with spaces
        arr[:] = " "

        for j in range(len(py_string)):
            # copy over each character
            arr[j] = py_string[j]

        return arr

    def _str_list_to_fort_char_array(self, strList, num_max_char):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takes a list of strings and returns the array.
        """

        arr = np.zeros((len(strList), num_max_char), dtype="str")
        arr[:] = " "
        for i, s in enumerate(strList):
            for j in range(len(s)):
                arr[i, j] = s[j]

        return arr

    def __fort_char_array_to_str(self, fort_string: str) -> str:
        # TODO: need a more general solution for |S<variable> type
        # SB: This should fix it but keep things commented out in case

        if fort_string.dtype == np.dtype("|S0"):
            # there are no characters in the sting to add
            return ""
        if fort_string.dtype == np.dtype("|S1"):
            py_string = b"".join(fort_string).decode().strip()
        elif fort_string.dtype == np.dtype("<U1"):
            py_string = "".join(fort_string).strip()
        elif fort_string.dtype.str.startswith("|S"):
            py_string = fort_string.decode().strip()
        # elif fort_string.dtype == np.dtype("|S4"):
        #     py_string = fort_string.decode().strip()
        # elif fort_string.dtype == np.dtype("|S16"):
        #     py_string = fort_string.decode().strip()
        # elif fort_string.dtype == np.dtype("|S40"):
        #     py_string = fort_string.decode().strip()
        # elif fort_string.dtype == np.dtype("|S80"):
        #     py_string = fort_string.decode().strip()
        # elif fort_string.dtype == np.dtype("|S120"):
        #     py_string = fort_string.decode().strip()
        else:
            raise TypeError(f"Unable to convert {fort_string} of type {fort_string.dtype} to string")

        return py_string

    def _fort_char_array_to_str_list(self, fortArray):
        """Undoes the _createFotranStringArray"""
        strList = []

        if fortArray.size == 1:
            # we must handle the 0-d array case sperately
            return self.__fort_char_array_to_str(fortArray[()])

        for ii in range(fortArray.size):
            py_string = self.__fort_char_array_to_str(fortArray[ii])
            if py_string != "":
                strList.append(py_string)

        return strList

    def _split_deriv_key(self, key):
        try:
            var, func = key.split("/")
        except Exception as e:
            print(key)
            raise e
        # remove leading 'd's
        var = var[1:]
        func = func[1:]
        return var, func

    def _get_deriv_key(self, var: str, func: str) -> str:
        return f"d{func}/d{var}"

    # ---------------------------
    # --- Derivative routines ---
    # ---------------------------

    # --- input ad seeds ---
    def get_variable_ad_seeds(self) -> Dict[str, float]:
        var_seeds = {}
        for con in self.con_var_to_fort_var:
            idx_con = self.conval_idx_dict[con]
            blk = "CASE_R" + self.ad_suffix
            var = "CONVAL" + self.ad_suffix
            slicer = (0, idx_con)

            fort_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)
            var_seeds[con] = copy.deepcopy(fort_arr)

        return var_seeds

    def set_variable_ad_seeds(self, con_seeds: Dict[str, Dict[str, float]], mode: str = "AD", scale=1.0) -> None:
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
                # # reverse lookup in the con_var_to_fort_var dict
                if con in self.con_surf_to_dindex:
                    val = self.get_control_deflection(con)
                else:
                    val = self.get_variable(con)

                val += con_seed_arr * scale

                if con in self.con_surf_to_dindex:
                    self.set_control_deflection(con, val)
                else:
                    self.set_variable(con, val)

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
            avl_key = self.ref_var_to_fort_var[ref_key]
            if len(avl_key) == 3:
                blk, var, slicer = avl_key
            else:
                blk, var = avl_key
                slicer = None

            if mode == "AD":
                blk += self.ad_suffix
                var += self.ad_suffix
                val = ref_seeds[ref_key] * scale
            elif mode == "FD":
                val = self.get_avl_fort_arr(blk, var, slicer=slicer)
                val += ref_seeds[ref_key] * scale

            self.set_avl_fort_arr(blk, var, val, slicer=slicer)

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
        for surf_key in self.unique_surface_names:
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
        var = "GAM_D_DIFF"

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
                cs_deriv_seeds[self._get_deriv_key(control, _var)] = val[()]

        return cs_deriv_seeds

    def set_consurf_derivs_ad_seeds(self, cs_deriv_seeds: Dict[str, float], scale=1.0):
        consurf_names = self.get_control_names()
        num_consurf = len(consurf_names)

        for deriv_func in cs_deriv_seeds:
            val_arr = np.zeros((num_consurf))

            var, cs_name = self._split_deriv_key(deriv_func)

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
            blk, var = self.case_stab_derivs_to_fort_var[func_key]

            blk += self.ad_suffix
            var += self.ad_suffix

            val = stab_deriv_seeds[func_key] * scale

            self.set_avl_fort_arr(blk, var, val)

    def get_body_axis_derivs_ad_seeds(self):
        deriv_data = {}

        for func_key, avl_vars in self.case_body_derivs_to_fort_var.items():
            deriv_data[func_key] = {}

            blk, var, slicer = avl_vars
            blk += self.ad_suffix
            var += self.ad_suffix
            val_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)
            deriv_data[func_key] = val_arr[()]

        return deriv_data

    def set_body_axis_derivs_ad_seeds(self, body_axis_deriv_seeds: Dict[str, Dict[str, float]], scale=1.0):
        for func_key in body_axis_deriv_seeds:
            blk, var, slicer = self.case_body_derivs_to_fort_var[func_key]

            blk += self.ad_suffix
            var += self.ad_suffix

            val = body_axis_deriv_seeds[func_key] * scale
            self.set_avl_fort_arr(blk, var, val, slicer=slicer)

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
        # use the size information to clear only the part of data that are used
        # zeroing more will cause more data to be allocated in physcial memory

        num_vor = self.get_mesh_size()
        gam = self.get_avl_fort_arr("VRTX_R", "GAM")
        num_vor_max = gam.size

        num_strips_max = self.get_avl_fort_arr("STRP_R", "AINC").size
        num_strips = self.get_avl_fort_arr("CASE_I", "NSTRIP")[()]

        num_surfs_max = self.get_avl_fort_arr("SURF_GEOM_R", "YDUPL").size
        num_surfs = self.get_avl_fort_arr("CASE_I", "NSURF")

        num_sec_max = self.get_avl_fort_arr("SURF_GEOM_R", "AINCS").shape[-1]
        num_sec = np.max(self.get_avl_fort_arr("SURF_GEOM_I", "NSEC"))

        num_airfoil_pts_max = self.get_avl_fort_arr("SURF_GEOM_R", "XLASEC").shape[-1]
        num_airfoil_pts = np.max(self.get_avl_fort_arr("SURF_GEOM_I", "NASEC"))

        # import psutil
        # process = psutil.Process()
        # mem_before = process.memory_info().rss
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

                            if dim_size == num_strips_max:
                                dim_size = num_strips

                            if dim_size == num_sec_max:
                                dim_size = num_sec

                            if dim_size == num_surfs_max:
                                dim_size = num_surfs

                            if dim_size == num_airfoil_pts_max:
                                dim_size = num_airfoil_pts

                            slices.append(slice(0, dim_size))
                        slicer = tuple(slices)
                        val[slicer] = 0.0

                        # mem_now = process.memory_info().rss
                        # del_mem = mem_now - mem_before
                        # mem_before = mem_now
                        # print(f'{att:10}:{_var:15}, {str(slicer):100} Memory usage: {del_mem/1024**2:.2f} MB')

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
    def _execute_jac_vec_prod_fwd(
        self,
        con_seeds: Optional[Dict[str, float]] = None,
        geom_seeds: Optional[Dict[str, Dict[str, any]]] = None,
        param_seeds: Optional[Dict[str, float]] = None,
        ref_seeds: Optional[Dict[str, float]] = None,
        gamma_seeds: Optional[np.ndarray] = None,
        gamma_d_seeds: Optional[np.ndarray] = None,
        gamma_u_seeds: Optional[np.ndarray] = None,
        mode: str = "AD",
        step: float = 1e-7,
    ) -> Tuple[Dict[str, float], np.ndarray, Dict[str, float], Dict[str, float], np.ndarray, np.ndarray]:
        """Get partial derivatives in forward mode. This routine is useful internally and when creating wrappers for things like OpenMDAO

        Args:
            con_seeds: Case constraint AD seeds
            geom_seeds: Geometric AD seeds in the same format as the geometric data
            param_seeds: Case parameter AD seeds
            ref_seeds: Reference condition AD seeds
            gamma_seeds: Circulation AD seeds
            gamma_d_seeds: dCirculation/d(Controls Deflection) AD seeds
            gamma_u_seeds:  dCirculation/d(flight condition) AD seeds
            mode: Either AD or FD. FD is mostly for testing
            step: Step size to use for the FD mode

        Returns:
            func_seeds: force coifficent AD seeds
            res_seeds:  residual AD seeds
            consurf_derivs_seeds: Control surface derivatives AD seeds
            stab_derivs_seeds: Stability derivatives AD seeds
            res_d_seeds: dResidual/d(Controls Deflection) AD seeds
            res_u_seeds: dResidual/d(flight condition) AD seeds
        """
        # TODO: add better name for gamma_d, it is too confusing

        # TODO: add error if data is not initailzed properly
        #   The easiest fix is to run an analysis or residual before hand

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
            self.set_variable_ad_seeds(con_seeds)
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
            body_axis_derivs_seeds = self.get_body_axis_derivs_ad_seeds()
            res_d_seeds = self.get_residual_d_ad_seeds()
            res_u_seeds = self.get_residual_u_ad_seeds()

            self.set_variable_ad_seeds(con_seeds, scale=0.0)
            self.set_geom_ad_seeds(geom_seeds, scale=0.0)
            self.set_gamma_ad_seeds(gamma_seeds, scale=0.0)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, scale=0.0)
            self.set_gamma_u_ad_seeds(gamma_u_seeds, scale=0.0)
            self.set_parameter_ad_seeds(param_seeds, scale=0.0)
            self.set_reference_ad_seeds(ref_seeds, scale=0.0)

            # TODO: remove??
            self.set_avl_fort_arr("VRTX_R_DIFF", "GAM_DIFF", gamma_seeds * 0.0, slicer=res_slice)

        if mode == "FD":
            self.set_variable_ad_seeds(con_seeds, mode="FD", scale=step)
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
            body_axis_deriv_petrub = self.get_body_axis_derivs()

            res_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES", slicer=res_slice))
            res_d_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_D", slicer=res_d_slice))
            res_u_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_U", slicer=res_u_slice))

            self.set_variable_ad_seeds(con_seeds, mode="FD", scale=-1 * step)
            self.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_ad_seeds(gamma_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_u_ad_seeds(gamma_u_seeds, mode="FD", scale=-1 * step)
            self.set_parameter_ad_seeds(param_seeds, mode="FD", scale=-1 * step)
            self.set_reference_ad_seeds(ref_seeds, mode="FD", scale=-1 * step)

            self.avl.update_surfaces()
            self.avl.get_res()
            self.avl.velsum()
            self.avl.aero()

            coef_data = self.get_total_forces()
            consurf_derivs = self.get_control_stab_derivs()
            stab_deriv = self.get_stab_derivs()
            body_axis_deriv = self.get_body_axis_derivs()

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
                stab_derivs_seeds[deriv_func] = (stab_deriv_petrub[deriv_func] - stab_deriv[deriv_func]) / step

            body_axis_derivs_seeds = {}
            for deriv_func in body_axis_deriv:
                body_axis_derivs_seeds[deriv_func] = (
                    body_axis_deriv_petrub[deriv_func] - body_axis_deriv[deriv_func]
                ) / step

            res_seeds = (res_peturbed - res) / step
            res_d_seeds = (res_d_peturbed - res_d) / step
            res_u_seeds = (res_u_peturbed - res_u) / step

        # TODO-clean: the way these arrays are returned is a bit of a mess
        return (
            func_seeds,
            res_seeds,
            consurf_derivs_seeds,
            stab_derivs_seeds,
            body_axis_derivs_seeds,
            res_d_seeds,
            res_u_seeds,
        )

    def _execute_jac_vec_prod_rev(
        self,
        func_seeds: Optional[Dict[str, float]] = None,
        res_seeds: Optional[np.ndarray] = None,
        consurf_derivs_seeds: Optional[Dict[str, float]] = None,
        stab_derivs_seeds: Optional[Dict[str, float]] = None,
        body_axis_derivs_seeds: Optional[Dict[str, float]] = None,
        res_d_seeds: Optional[np.ndarray] = None,
        res_u_seeds: Optional[np.ndarray] = None,
        print_timings: Optional[bool] = False,
    ) -> Tuple[
        Dict[str, float],
        Dict[str, Dict[str, any]],
        np.ndarray,
        np.ndarray,
        np.ndarray,
        Dict[str, float],
        Dict[str, float],
    ]:
        """Get partial derivatives in reverse mode. This routine is useful internally and when creating wrappers for things like OpenMDAO

        Args:
            func_seeds: force coefficient AD seeds
            res_seeds:  residual AD seeds
            consurf_derivs_seeds: Control surface derivatives AD seeds
            stab_derivs_seeds: Stability derivatives AD seeds
            body_axis_derivs_seeds: body axis derivatives AD seeds
            res_d_seeds: dResidual/d(Controls Deflection) AD seeds
            res_u_seeds: dResidual/d(flight condition) AD seeds
            print_timings: flag to show timing data

        Returns:
            con_seeds: Case constraint AD seeds
            geom_seeds: Geometric AD seeds in the same format as the geometric data
            gamma_seeds: Circulation AD seeds
            gamma_d_seeds: dCirculation/d(Controls Deflection) AD seeds
            gamma_u_seeds:  dCirculation/d(flight condition) AD seeds
            param_seeds: Case parameter AD seeds
            ref_seeds: Reference condition AD seeds
        """
        # extract derivatives seeds and set the output dict of functions

        if print_timings:
            time_start = time.time()

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

        if body_axis_derivs_seeds is None:
            body_axis_derivs_seeds = {}

        # set derivative seeds
        # self.clear_ad_seeds()
        time_last = time.time()
        self.set_function_ad_seeds(func_seeds)
        self.set_residual_ad_seeds(res_seeds)
        self.set_residual_d_ad_seeds(res_d_seeds)
        self.set_residual_u_ad_seeds(res_u_seeds)
        self.set_consurf_derivs_ad_seeds(consurf_derivs_seeds)
        self.set_stab_derivs_ad_seeds(stab_derivs_seeds)
        self.set_body_axis_derivs_ad_seeds(body_axis_derivs_seeds)

        if print_timings:
            print(f"    Time to set seeds: {time.time() - time_last}")
            time_last = time.time()

        # propogate the seeds through without resolveing
        self.avl.aero_b()
        if print_timings:
            print(f"    Time to propogate seeds:aero_b: {time.time() - time_last}")
            time_last = time.time()
        self.avl.velsum_b()
        if print_timings:
            print(f"    Time to propogate seeds:velsum_b: {time.time() - time_last}")
            time_last = time.time()
        self.avl.get_res_b()
        if print_timings:
            print(f"    Time to propogate seeds:get_res_b: {time.time() - time_last}")
            time_last = time.time()
        self.avl.update_surfaces_b()
        if print_timings:
            print(f"    Time to propogate seeds:update_surfaces_b: {time.time() - time_last}")
            time_last = time.time()

        # extract derivatives seeds and set the output dict of functions
        con_seeds = self.get_variable_ad_seeds()
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
        self.set_body_axis_derivs_ad_seeds(body_axis_derivs_seeds, scale=0.0)
        if print_timings:
            print(f"    Time to clear seeds: {time.time() - time_last}")
            time_last = time.time()

        if print_timings:
            print(f"   Total Time: {time.time() - time_start}")

        return con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds, gamma_u_seeds, param_seeds, ref_seeds

    def execute_run_sensitivities(
        self,
        funcs: List[str],
        stab_derivs: Optional[List[str]] = None,
        body_axis_derivs: Optional[List[str]] = None,
        consurf_derivs: Optional[List[str]] = None,
        print_timings: Optional[bool] = False,
    ) -> Dict[str, Dict[str, float]]:
        """Run the sensitivities of the input functionals in adjoint mode

        Args:
            funcs: force coefficients to compute the sensitivities with respect to
            stab_derivs: stability derivatives to compute the sensitivities with respect to
            body_axis_derivs: body axis derivatives to compute the sensitivities with respect to
            consurf_derivs: control surface derivates to compute the sensitivities with respect to
            print_timings: flag to print timing information

        Returns:
            sens: a nested dictionary of sensitivities. The first key is the function and the next keys are for the design variables.
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
            _, _, pfpU, _, _, _, _ = self._execute_jac_vec_prod_rev(func_seeds={func: 1.0})
            if print_timings:
                print(f"Time to get RHS: {time.time() - time_last}")
                time_last = time.time()

            # self.clear_ad_seeds()
            # u solver adjoint equation with RHS
            self.set_gamma_ad_seeds(-1 * pfpU)
            solve_gamma_u_adj = False
            solve_gamma_d_adj = False
            self.avl.solve_adjoint(solve_gamma_u_adj, solve_gamma_d_adj)
            if print_timings:
                print(f"Time to solve adjoint: {time.time() - time_last}")
                time_last = time.time()
            # get the resulting adjoint vector (dfunc/dRes) from fortran
            dfdR = self.get_residual_ad_seeds()
            # self.clear_ad_seeds()
            con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self._execute_jac_vec_prod_rev(
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
                _, _, pfpU, pf_pU_d, _, _, _ = self._execute_jac_vec_prod_rev(consurf_derivs_seeds={func_key: 1.0})
                if print_timings:
                    print(f"Time to get RHS: {time.time() - time_last}")
                    time_last = time.time()

                # self.clear_ad_seeds()
                # u solver adjoint equation with RHS
                self.set_gamma_ad_seeds(-1 * pfpU)
                self.set_gamma_d_ad_seeds(-1 * pf_pU_d)
                solve_gamma_u_adj = False
                solve_gamma_d_adj = True
                self.avl.solve_adjoint(solve_gamma_u_adj, solve_gamma_d_adj)
                if print_timings:
                    print(f"Time to solve adjoint: {time.time() - time_last}")
                    time_last = time.time()

                # get the resulting adjoint vector (dfunc/dRes) from fortran
                dfdR = self.get_residual_ad_seeds()
                dfdR_d = self.get_residual_d_ad_seeds()
                # self.clear_ad_seeds()
                con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self._execute_jac_vec_prod_rev(
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

            for func_key in stab_derivs:
                if func_key not in sens:
                    sens[func_key] = {}

                # get the RHS of the adjoint equation (pFpU)
                # TODO: remove seeds if it doesn't effect accuracy
                _, _, pfpU, _, pf_pU_u, _, _ = self._execute_jac_vec_prod_rev(stab_derivs_seeds={func_key: 1.0})
                if print_timings:
                    print(f"Time to get RHS: {time.time() - time_last}")
                    time_last = time.time()

                # self.clear_ad_seeds()
                # u solver adjoint equation with RHS
                self.set_gamma_ad_seeds(-1 * pfpU)
                self.set_gamma_u_ad_seeds(-1 * pf_pU_u)
                solve_gamma_u_adj = True
                solve_gamma_d_adj = False
                self.avl.solve_adjoint(solve_gamma_u_adj, solve_gamma_d_adj)
                if print_timings:
                    print(f"Time to solve adjoint: {time.time() - time_last}")
                    time_last = time.time()

                # get the resulting adjoint vector (dfunc/dRes) from fortran
                dfdR = self.get_residual_ad_seeds()
                dfdR_u = self.get_residual_u_ad_seeds()
                # self.clear_ad_seeds()
                con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self._execute_jac_vec_prod_rev(
                    stab_derivs_seeds={func_key: 1.0}, res_seeds=dfdR, res_u_seeds=dfdR_u
                )

                if print_timings:
                    print(f"Time to combine : {time.time() - time_last}")
                    time_last = time.time()

                sens[func_key].update(con_seeds)
                sens[func_key].update(geom_seeds)
                sens[func_key].update(param_seeds)
                sens[func_key].update(ref_seeds)
                # sd_deriv_seeds[func_key] = 0.0

        if body_axis_derivs is not None:
            if print_timings:
                print("Running stab derivs")
                time_last = time.time()

            for func_key in body_axis_derivs:
                if func_key not in sens:
                    sens[func_key] = {}

                # get the RHS of the adjoint equation (pFpU)
                # TODO: remove seeds if it doesn't effect accuracy
                _, _, pfpU, _, pf_pU_u, _, _ = self._execute_jac_vec_prod_rev(body_axis_derivs_seeds={func_key: 1.0})
                if print_timings:
                    print(f"Time to get RHS: {time.time() - time_last}")
                    time_last = time.time()

                # self.clear_ad_seeds()
                # u solver adjoint equation with RHS
                self.set_gamma_ad_seeds(-1 * pfpU)
                self.set_gamma_u_ad_seeds(-1 * pf_pU_u)
                solve_gamma_u_adj = True
                solve_gamma_d_adj = False
                self.avl.solve_adjoint(solve_gamma_u_adj, solve_gamma_d_adj)
                if print_timings:
                    print(f"Time to solve adjoint: {time.time() - time_last}")
                    time_last = time.time()

                # get the resulting adjoint vector (dfunc/dRes) from fortran
                dfdR = self.get_residual_ad_seeds()
                dfdR_u = self.get_residual_u_ad_seeds()
                # self.clear_ad_seeds()
                con_seeds, geom_seeds, _, _, _, param_seeds, ref_seeds = self._execute_jac_vec_prod_rev(
                    body_axis_derivs_seeds={func_key: 1.0}, res_seeds=dfdR, res_u_seeds=dfdR_u
                )

                if print_timings:
                    print(f"Time to combine : {time.time() - time_last}")
                    time_last = time.time()

                sens[func_key].update(con_seeds)
                sens[func_key].update(geom_seeds)
                sens[func_key].update(param_seeds)
                sens[func_key].update(ref_seeds)

        return sens

    # --- ploting and vizulaization ---
    def add_mesh_plot(
        self,
        axis,
        xaxis: str = "x",
        yaxis: str = "y",
        color: str = "black",
        mesh_style="--",
        mesh_linewidth=0.3,
        show_mesh: bool = True,
    ):
        """Adds a plot of the aircraft mesh to the axis

        Args:
            axis: axis to add the plot to
            xaxis: what variable should be plotted on the x axis. Options are ['x', 'y', 'z']
            yaxis: what variable should be plotted on the y-axis. Options are ['x', 'y', 'z']
            color: what color should the mesh be
            mesh_style: line style of the interior mesh, e.g. '-' or '--'
            mesh_linewidth: width of the interior mesh, 1.0 will match the surface outline
            show_mesh: flag to show the interior mesh of the geometry
        """
        mesh_size = self.get_mesh_size()
        num_control_surfs = self.get_num_control_surfs()
        num_strips = self.get_num_strips()
        num_surfs = self.get_num_surfaces()

        # get the mesh points for ploting
        mesh_slice = (slice(0, mesh_size),)
        strip_slice = (slice(0, num_strips),)
        surf_slice = (slice(0, num_surfs),)

        rv1 = self.get_avl_fort_arr("VRTX_R", "RV1", slicer=mesh_slice)  # Vortex Left points
        rv2 = self.get_avl_fort_arr("VRTX_R", "RV2", slicer=mesh_slice)  # Vortex Right points
        rle1 = self.get_avl_fort_arr("STRP_R", "RLE1", slicer=strip_slice)  # Strip left end LE point
        rle2 = self.get_avl_fort_arr("STRP_R", "RLE2", slicer=strip_slice)  # Strip right end LE point
        chord1 = self.get_avl_fort_arr("STRP_R", "CHORD1", slicer=strip_slice)  # Left strip chord
        chord2 = self.get_avl_fort_arr("STRP_R", "CHORD2", slicer=strip_slice)  # Right strip chord
        jfrst = self.get_avl_fort_arr("SURF_I", "JFRST", slicer=surf_slice)  # Index of first strip in surface

        ijfrst = self.get_avl_fort_arr("STRP_I", "IJFRST", slicer=strip_slice)  # Index of first element in strip
        nvstrp = self.get_avl_fort_arr("STRP_I", "NVSTRP", slicer=strip_slice)  # Number of elements in strip

        nj = self.get_avl_fort_arr("SURF_I", "NJ", slicer=surf_slice)  # Number of elements along span in surface
        imags = self.get_avl_fort_arr("SURF_I", "IMAGS")  # Is surface YDUPL one?

        for idx_surf in range(num_surfs):
            # get the range of the elements that belong to this surfaces
            strip_st = jfrst[idx_surf] - 1
            strip_end = strip_st + nj[idx_surf]

            # inboard and outboard of outline
            # get surfaces that have not been duplicated
            if imags[idx_surf] > 0:
                j1 = strip_st
                jn = strip_end - 1
                dj = 1
            else:
                # this surface is a duplicate
                j1 = strip_end - 1
                jn = strip_st
                dj = -1

            pts = {
                "x": [rle1[j1, 0], rle1[j1, 0] + chord1[j1]],
                "y": [rle1[j1, 1], rle1[j1, 1]],
                "z": [rle1[j1, 2], rle1[j1, 2]],
            }
            # # chord-wise grid
            axis.plot(pts[xaxis], pts[yaxis], color=color)

            pts = {
                "x": np.array([rle2[jn, 0], rle2[jn, 0] + chord2[jn]]),
                "y": np.array([rle2[jn, 1], rle2[jn, 1]]),
                "z": np.array([rle2[jn, 2], rle2[jn, 2]]),
            }

            # # chord-wise grid
            axis.plot(pts[xaxis], pts[yaxis], color=color)

            # # --- outline of surface ---
            # front
            pts = {
                "x": np.append(rle1[j1:jn:dj, 0], rle2[jn, 0]),
                "y": np.append(rle1[j1:jn:dj, 1], rle2[jn, 1]),
                "z": np.append(rle1[j1:jn:dj, 2], rle2[jn, 2]),
            }
            axis.plot(pts[xaxis], pts[yaxis], "-", color=color)

            # aft

            pts = {
                "x": np.append(rle1[j1:jn:dj, 0] + chord1[j1:jn:dj], rle2[jn, 0] + chord2[jn]),
                "y": np.append(rle1[j1:jn:dj, 1], rle2[jn, 1]),
                "z": np.append(rle1[j1:jn:dj, 2], rle2[jn, 2]),
            }
            axis.plot(pts[xaxis], pts[yaxis], "-", color=color)

            if show_mesh:
                for idx_strip in range(strip_st, strip_end):
                    if ((imags[idx_surf] > 0) and (idx_strip != strip_st)) or (
                        (imags[idx_surf] < 0) and (idx_strip != strip_end)
                    ):
                        pts = {
                            "x": [rle1[idx_strip, 0], rle1[idx_strip, 0] + chord1[idx_strip]],
                            "y": [rle1[idx_strip, 1], rle1[idx_strip, 1]],
                            "z": [rle1[idx_strip, 2], rle1[idx_strip, 2]],
                        }

                        # # chord-wise grid
                        axis.plot(pts[xaxis], pts[yaxis], mesh_style, color=color, alpha=0.7, linewidth=mesh_linewidth)

                    vor_st = ijfrst[idx_strip] - 1
                    vor_end = vor_st + nvstrp[idx_strip]

                    # spanwise grid
                    for idx_vor in range(vor_st, vor_end):
                        pts = {
                            "x": [rv1[idx_vor, 0], rv2[idx_vor, 0]],
                            "y": [rv1[idx_vor, 1], rv2[idx_vor, 1]],
                            "z": [rv1[idx_vor, 2], rv2[idx_vor, 2]],
                        }
                        axis.plot(pts[xaxis], pts[yaxis], mesh_style, color=color, alpha=0.7, linewidth=mesh_linewidth)

    def add_body_plot(
        self,
        axis,
        xaxis: str = "x",
        yaxis: str = "y",
        color: str = "magenta",
    ):
        """Adds a plot of the body geometry to the axis

        Args:
            axis: axis to add the plot to
            xaxis: what variable should be plotted on the x axis. Options are ['x', 'y', 'z']
            yaxis: what variable should be plotted on the y-axis. Options are ['x', 'y', 'z']
            color: what color should the body be plotted in
        """
        import numpy as np

        # Check if any bodies exist
        num_bodies = self.get_avl_fort_arr("CASE_I", "NBODY")
        if num_bodies == 0:
            return

        # Get body data
        nl = self.get_avl_fort_arr("BODY_I", "NL", slicer=slice(0, num_bodies))
        lfrst = self.get_avl_fort_arr("BODY_I", "LFRST", slicer=slice(0, num_bodies))
        print(lfrst)

        # Calculate total number of body nodes needed
        # We need lfrst[i] + nl[i] to include all nodes we might access
        max_node_idx = max(lfrst[i] + nl[i] for i in range(num_bodies))

        # Get body line coordinates and radii
        rl = self.get_avl_fort_arr("VRTX_R", "RL", slicer=(slice(None), slice(0, max_node_idx)))
        radl = self.get_avl_fort_arr("VRTX_R", "RADL", slicer=slice(0, max_node_idx))

        # Number of points around each cross-section circle
        kk = 51

        # Plot each body
        for n in range(num_bodies):
            # First and last node indices for this body (Fortran uses 1-based indexing)
            l1 = lfrst[n] - 1  # Convert to 0-based
            ln = l1 + nl[n]

            # Plot centerline
            centerline_pts = {
                "x": rl[l1:ln, 0],
                "y": rl[l1:ln, 1],
                "z": rl[l1:ln, 2],
            }
            axis.plot(centerline_pts[xaxis], centerline_pts[yaxis], "-", color=color, linewidth=1.5)

            # Plot cross-section circles at each node
            for l in range(l1, ln):
                # Calculate tangent vector from neighboring nodes
                lm = max(l - 1, l1)
                lp = min(l + 1, ln - 1)

                dxl = rl[lp,0] - rl[lm, 0]
                dyl = rl[lp,1] - rl[lm, 1]
                dzl = rl[lp,2] - rl[lm, 2]
                dsl = np.sqrt(dxl**2 + dyl**2 + dzl**2)

                if dsl != 0.0:
                    dxl /= dsl
                    dyl /= dsl
                    dzl /= dsl

                # Calculate perpendicular unit vectors
                # Horizontal vector (in x-y plane, perpendicular to tangent)
                uhx = -dyl
                uhy = dxl
                uhz = 0.0

                # Vertical vector (perpendicular to both tangent and horizontal)
                uvx = dyl * uhz - dzl * uhy
                uvy = dzl * uhx - dxl * uhz
                uvz = dxl * uhy - dyl * uhx

                # Generate circle points
                circle_x = []
                circle_y = []
                circle_z = []

                for k in range(kk + 1):  # +1 to close the circle
                    theta = 2.0 * np.pi * k / kk

                    # Components of the circle point
                    hx = uhx * np.cos(theta)
                    hy = uhy * np.cos(theta)
                    hz = uhz * np.cos(theta)

                    vx = uvx * np.sin(theta)
                    vy = uvy * np.sin(theta)
                    vz = uvz * np.sin(theta)

                    # Circle point at radius RADL
                    circle_x.append(rl[l, 0] + (hx + vx) * radl[l])
                    circle_y.append(rl[l, 1] + (hy + vy) * radl[l])
                    circle_z.append(rl[l, 2] + (hz + vz) * radl[l])

                # Plot the circle
                circle_pts = {
                    "x": np.array(circle_x),
                    "y": np.array(circle_y),
                    "z": np.array(circle_z),
                }
                axis.plot(circle_pts[xaxis], circle_pts[yaxis], "-", color=color, linewidth=0.5)

    def plot_geom(self, axes=None, body_color="magenta"):
        """Generate a matplotlib plot of geometry

        Args:
            axes: Matplotlib axis object to add the plots too. If none are given, the axes will be generated.
            body_color: Color to use for plotting body geometry (default: 'magenta')
        """

        if axes == None:
            import matplotlib.pyplot as plt

            ax1 = plt.subplot(2, 1, 1)
            ax2 = plt.subplot(2, 1, 2)

            ax2.set_ylabel("Z", rotation=0)
            ax2.set_xlabel("Y")
            ax1.set_ylabel("X", rotation=0)
            ax1.set_aspect("equal")
            ax2.set_aspect("equal")

        else:
            ax1, ax2 = axes

        self.add_mesh_plot(ax1, xaxis="y", yaxis="x")
        self.add_body_plot(ax1, xaxis="y", yaxis="x", color=body_color)

        self.add_mesh_plot(ax2, xaxis="y", yaxis="z")
        self.add_body_plot(ax2, xaxis="y", yaxis="z", color=body_color)

        if axes == None:
            # assume that if we don't provide axes that we want to see the plot
            plt.axis("equal")
            plt.show()

    def get_cp_data(self) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """Gets the current surface mesh and cp distribution

        Returns:
            xyz_list: list of surface mesh points
            cp_list: list of cp points
        """
        nvort = self.get_mesh_size()
        mesh_slicer = (slice(0, nvort),)

        num_surfs = self.get_num_surfaces()
        surf_slice = slice(0, num_surfs)

        nj = self.get_avl_fort_arr("SURF_I", "NJ", slicer=surf_slice)
        nk = self.get_avl_fort_arr("SURF_I", "NK", slicer=surf_slice)

        xyz_list = []
        cp_list = []

        for idx_surf in range(num_surfs):
            nChords = nk[idx_surf]
            nStrips = nj[idx_surf]
            nPts = (nStrips + 1) * (nChords * 2 + 1)
            nCCPts = (nStrips) * (nChords * 2)

            vtx_slicer = (idx_surf, slice(0, nPts), slice(None))
            cp_slicer = (idx_surf, slice(0, nCCPts))

            xyz = self.get_avl_fort_arr("VRTX_S", "XYZSURF", slicer=vtx_slicer)

            xyz = xyz.reshape((nStrips + 1, nChords * 2 + 1, 3))
            xyz_list.append(xyz)

            cp = self.get_avl_fort_arr("VRTX_S", "CPSURF", slicer=cp_slicer)
            cp = cp.reshape((nStrips, nChords * 2))
            cp_list.append(cp)

        return xyz_list, cp_list

    def plot_cp(self):
        """Create a matplotlib plot of the surface and cp distribution"""
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
        norm = plt.Normalize(vmin=-1 * cp_amax, vmax=cp_amax)
        m = cm.ScalarMappable(cmap=cm.bwr)
        m.set_clim(vmin=cp_min, vmax=cp_max)

        # do the actual ploting of each surface
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        plt.subplots_adjust(left=0.025, right=0.925, top=0.925, bottom=0.025)

        for idx_surf in range(num_surfs):
            xyz_surf = xyz[idx_surf]

            face_color = m.to_rgba(cp[idx_surf])

            ax.plot_surface(xyz_surf[:, :, 0], xyz_surf[:, :, 1], xyz_surf[:, :, 2], facecolors=face_color)

        plt.axis("off")

        plt.grid(b=None)
        colorbar = fig.colorbar(m, ax=ax)
        colorbar.set_label("Cp", rotation=0, labelpad=20)
        # Set an equal aspect ratio
        ax.set_aspect("equal")

        plt.show()

    def write_tecplot(self, file_name: str, solution_time: float = None):
        """Write a tecplot file of the current surface and Cp distribution

        Args:
            file_name: Name of the output file
            solution_time: Add a solution time to the output. This is useful for flipping through data in tecplot, but breaks Paraview.
        """
        if solution_time is not None:
            add_time = True
        else:
            solution_time = 0.0
            add_time = False

        self.avl.cpoml(False)
        self.avl.write_tecplot(file_name + ".dat", add_time, solution_time)
