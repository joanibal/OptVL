import numpy as np


def get_ovl_input_name(avl_input, ovl):
    avl_to_ovl_vars = {
        "A lpha": "alpha",
        "alpha": "alpha",
        "B eta": "beta",
        "beta": "beta",
        "R oll  rate": "roll rate",
        "pb/2V": "roll rate",
        "P itch rate": "pitch rate",
        "qc/2V": "pitch rate",
        "Y aw   rate": "yaw rate",
        "rb/2V": "yaw rate",
        "CL": "CL",
        "CY": "CY",
        "Cl roll mom": "Cl",
        "Cm pitchmom": "Cm",
        "Cn yaw  mom": "Cn",
    }
    if avl_input in avl_to_ovl_vars:
        return avl_to_ovl_vars[avl_input]
    else:
        for consurf in ovl.get_control_names():
            if consurf in avl_input:
                return consurf

        raise ValueError(f"avl input {avl_input} could not be mapped to OptVL input")


def get_avl_output_name(ovl_output, ovl):
    ovl_to_avl_vars = {
        "Sref": "Sref",
        "Cref": "Cref",
        "Bref": "Bref",
        "Xref": "Xref",
        "Yref": "Yref",
        "Zref": "Zref",
        "alpha": "Alpha",
        "roll rate": "pb/2V",
        #   "p'b/2V",
        "beta": "Beta",
        "pitch rate": "qc/2V",
        "mach": "Mach",
        "yaw rate": "rb/2V",
        #   "r'b/2V",
        "CX body axis": "CXtot",
        "CY body axis": "CYtot",
        "CZ body axis": "CZtot",
        "CX": "CXtot",
        "CY": "CYtot",
        "CZ": "CZtot",
        "CL": "CLtot",
        "CD": "CDtot",
        "CDv": "CDvis",
        "CDi": "CDind",
        "Cl": "Cltot",
        "Cl'": "Cl'tot",
        "Cm": "Cmtot",
        "Cn": "Cntot",
        "Cn'": "Cn'tot",
        "CLff": "CLff",
        "CDff": "CDff",
        "CYff": "CYff",
        "e": "e",
    }
    if ovl_output in ovl_to_avl_vars:
        return ovl_to_avl_vars[ovl_output]
    else:
        for consurf in ovl.get_control_names():
            if consurf == ovl_output:
                return consurf

        raise ValueError(f"ovl output {ovl_output} could not be mapped to AVL Output")


def set_inputs(ovl, ref_data):
    input_dict = ref_data["inputs"]

    for var_dict in input_dict:
        avl_var = var_dict["variable"].strip()
        avl_con = var_dict["constraint"].strip()
        val = var_dict["value"]

        # we need to translate the AVL data into OptVL inputs
        var = get_ovl_input_name(avl_var, ovl)
        con = get_ovl_input_name(avl_con, ovl)

        ovl.set_constraint(var, con, val)


# TODO: change the partial derivs test to use this to make it easier to
# debug via printing
def check_vals(
    val,
    ref_val,
    msg,
    rtol=1e-15,
    atol=0,
    printing=False,
    raise_error=True,
    colored_print=True,
    label_width=10,
    val_fmt="13.5e",
    err_fmt="11.3e",
):
    if printing:
        abs_err = val - ref_val
        if ref_val != 0.0:
            rel_err = abs_err / ref_val
        else:
            rel_err = np.nan

        # if the tolerances are violated print it in red
        test_failed = np.abs(abs_err) > (atol + rtol * np.abs(ref_val))

        # if the tolerances are violated print it in red
        if colored_print and test_failed:
            start = "\033[91m"
            end = "\033[0m"
        elif test_failed:
            start = "X "
            end = ""
        else:
            if colored_print:
                start = ""
                end = ""
            else:
                start = "  "
                end = ""

        print(
            f"{start}{msg:{label_width}} {val:{val_fmt}} {ref_val:{val_fmt}} | {rel_err:{err_fmt}} {abs_err:{err_fmt}}{end}"
        )

    if raise_error:
        np.testing.assert_allclose(
            val,
            ref_val,
            rtol=rtol,
            atol=atol,
            err_msg=msg,
        )


def run_comparison(ovl, ref_data_cases, **kwargs):
    for ref_data in ref_data_cases:
        set_inputs(ovl, ref_data)

        # This is the method that acutally runs the analysis
        ovl.execute_run()

        mesh_data = ovl.get_mesh_data()
        for key in ref_data["outputs"]["mesh_data"]:
            mesh_avl_val = ref_data["outputs"]["mesh_data"][key]
            check_vals(mesh_data[key], mesh_avl_val, f"meshing:{key}", **kwargs)

        reference_data = ovl.get_reference_data()
        for key in reference_data:
            ref_avl_val = ref_data["outputs"]["reference_data"][key]
            check_vals(reference_data[key], ref_avl_val, f"ref:{key}", **kwargs)

        force_data = ovl.get_total_forces()
        for key in force_data:
            avl_key = get_avl_output_name(key, ovl)
            avl_val = ref_data["outputs"]["total_forces"][avl_key]
            check_vals(force_data[key], avl_val, key, **kwargs)

        surf_avl_keys = {
            'length':"Length",
            'area': "Area",
            # 'area': "Ssurf",
            'volume': "Vol",
            "average chord":"Cave",
            "CL surf" : "cl",
            "CD surf" : "cd",
            "CDv surf" : "cdv",
        }

        surface_data = ovl.get_surface_forces()
        for idx_surf, surf in enumerate(surface_data):
            for key in surface_data[surf]:
                
                if key in surf_avl_keys:
                    avl_key = surf_avl_keys[key]
                elif key in ["CX", "CY", "CZ", "CMLE_LSTRP surf"]:
                    # avl output does not provide this data
                    continue
                else:
                    avl_key = key

                avl_val = ref_data["outputs"]["surface_forces"][surf][avl_key]
                check_vals(surface_data[surf][key], avl_val, key, **kwargs)

        strip_avl_keys = {
            'X LE':"Xle",
            'Y LE':"Yle",
            'Z LE':"Zle",
            'chord':"Chord",
            'area': "Area",
            'spanloading': "c_cl",
            "downwash":"ai",
            "CL perp" : "cl_perp",
            "CL strip" : "cl",
            "CD strip" : "cd",
            "CDv" : "cdv",
            "Cm c/4" : "cm_c/4",
            "Cm LE" : "cm_LE",
            "CP x/c" : "C.P.x/c",
        }

        strip_data = ovl.get_strip_forces()
        
        
        for idx_surf, surf in enumerate(strip_data):
            for key in strip_data[surf]:
            
                
                if key in strip_avl_keys:
                    avl_key = strip_avl_keys[key]
                elif key in ["S LE", "width", "twist", "CL", "CD", "CX", "CY", "CZ", "Cl", "Cm", "Cn", "CF strip", "Cm strip", "CP x/c", "lift dist", "drag dist", "roll dist", "yaw dist"]:
                    # avl output does not provide this data
                    continue
                else:
                    avl_key = key

                
                avl_val = ref_data["outputs"]["strip_forces"][surf][avl_key]
                check_vals(strip_data[surf][key], avl_val, key, **kwargs)

        
        body_avl_keys = {
            'length':"Length",
            'surface area': "Asurf",
            'volume': "Vol",
        }
        body_force_data = ovl.get_body_forces()
        for idx_body, body in enumerate(body_force_data):
            for key in body_force_data[body]:
                if key in body_avl_keys:
                    avl_key = body_avl_keys[key]
                else:
                    avl_key = key
                avl_val = ref_data["outputs"]["body_forces"][str(idx_body+1)][avl_key]
                check_vals(body_force_data[body][key], avl_val, key, **kwargs)

        control_def = ovl.get_control_deflections()
        for key in control_def:
            def_avl_val = ref_data["outputs"]["controls"][key]
            check_vals(control_def[key], def_avl_val, key, **kwargs)

        stab_derivs = ovl.get_stab_derivs()
        for key in stab_derivs:
            # we must translate the key into the AVL key
            if key == "lateral parameter":
                continue
            if key == "static margin":
                continue
            elif key == "neutral point":
                avl_key = "Xnp"
            elif key == "spiral parameter":
                avl_key = key
                # this behavior is modified relative to AVL
                if np.abs(stab_derivs["dCl'/dr'"] * stab_derivs["dCn'/dbeta"]) <= 0.0001:
                    continue
            else:
                func, var = key[1:].split("/d")
                avl_key = f"{func[:2]}{var[0]}"

            deriv_avl_val = ref_data["outputs"]["stability_derivatives"][avl_key]
            
            check_vals(stab_derivs[key], deriv_avl_val, f"sd:{key}", **kwargs)

        body_axi_derivs = ovl.get_body_axis_derivs()
        for key in body_axi_derivs:
            # we must translate the key into the AVL key
            func, var = key[1:].split("/d")
            avl_key = f"{func[:2]}{var[0]}"

            deriv_avl_val = ref_data["outputs"]["body_axis_derivatives"][avl_key]
            
            check_vals(body_axi_derivs[key], deriv_avl_val, f"bd:{key}", **kwargs)

        con_stab_derivs = ovl.get_control_stab_derivs()
        for key in con_stab_derivs:
            func, var = key[1:].split("/d")

            avl_key = f"{func[:2]}d{var}"

            if func in ["CL", "CD", "Cl'", "Cn'"]:
                axis_key = "stability_derivatives"
            else:
                axis_key = "body_axis_derivatives"

            avl_val = ref_data["outputs"][axis_key][avl_key]
            check_vals(con_stab_derivs[key], avl_val, f"cd:{key}", **kwargs)
