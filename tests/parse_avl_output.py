import re
import json
import argparse

# --- Regex patterns ---
line_re = re.compile(
    r"^\s*([A-Za-z0-9'\-/ ]+?)\s*->\s*([A-Za-z0-9'\-/ ]+?)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[Ee][-+]?\d+)?)",
    re.M
)

number_re = re.compile(r"([-+]?[0-9]*\.?[0-9]+(?:[Ee][-+]?\d+)?)")

# --- Helper to safely parse floats ---
def parse_floats(line):
    return [float(x) for x in number_re.findall(line)]

# --- Main parser ---
def parse_avl_output(run_output_files):
    cases = []
    case_idx = 0
    current_case = None
    
    with open(run_output_files) as fid:
        file_text = fid.read()
    
        lines = iter(file_text.splitlines())
        for line in lines:
            line = line.strip()
            
            
            # save the input data 
            if "variable          constraint" in line: 
                input_data = []
                
                next(lines) # skip header
                while True:
                    l = next(lines)
                    
                    if '------------' in l:
                        # we are done 
                        break
                    
                    var, con = l.split('->')
                    con, con_val = con.split('=')
                    
                    input_data.append({
                        'variable':var.strip(),
                        'constraint':con.strip(),
                        'value': float(con_val)
                    })
            

            # --- Detect start of output block ---
            if "TOT" == line:
                
                if current_case:
                    cases.append(current_case)
                case_idx += 1
                current_case = {
                    "case_number": case_idx,
                    "inputs": input_data,
                    "outputs": {
                        "mesh_data": {},
                        "reference_data": {},
                        "variables": {},
                        "total_forces": {},
                        "surface_forces": {},
                        "strip_forces": {},
                        "body_forces": {},
                        "controls": {},
                        "stability_derivatives": {},
                        "body_axis_derivatives": {}
                    }
                }
                print(f'adding case, {case_idx}')
                
                # skip the next 3 lines
                for _ in range(3):
                    next(lines)
                # continue

                # --- Mesh data ---
                for _ in range(3):
                    val, var = next(lines).split('|')
                    # remove '#'
                    var = var.replace('#', '')
                    current_case["outputs"]["mesh_data"][var.strip()] = int(val)

                # --- Reference data ---
                # Sref Cref Bref
                vals, keys = next(lines).split("|")
                vals = parse_floats(vals)
                keys = keys.split(',')
                
                for k, v in zip(keys, vals):
                    current_case["outputs"]["reference_data"][k.strip()] = v
                
                # Xref Yref Zref
                vals, keys = next(lines).split("|")
                vals = parse_floats(vals)
                keys = keys.split(',')
                current_case["outputs"]["reference_data"]['XYZref'] = vals
                

                # skip 2
                for _ in range(2):
                    next(lines) 
                    
                # variables
                for _ in range(3):
                    vals, keys = next(lines).split("|")
                    vals = parse_floats(vals)
                    keys = keys.split(',')
                    
                    for k, v in zip(keys, vals):
                        current_case["outputs"]["variables"][k.strip()] = v
                    
                # --- Coefficients (total forces) ---
                for _ in range(7):
                    vals, keys = next(lines).split("|")
                    vals = parse_floats(vals)
                    keys = keys.split(':')[-1]
                    keys = keys.split(',')
                    
                    for k, v in zip(keys, vals):
                        current_case["outputs"]["total_forces"][k.strip()] = v

                # --- Controls ---
                next(lines) # skip the header
                ncontrols = int(next(lines))
                for _ in range(ncontrols):
                    l = next(lines).strip()
                    val, name  = l.split("  ")
                    current_case["outputs"]["controls"][name.strip()] = float(val)
            
            elif "s>  BODY" in line:
                # skip version line
                # and reference info
                for _ in range(4):
                    next(lines)                # skip a
                
                nbodys, _ = next(lines).split("|")
                nbodys = int(nbodys)
                
                # skip next two lines
                for idx_body in range(1,nbodys+1):
                    for _ in range(1):
                        next(lines)
                    
                    body_key = next(lines).strip()
                    # For some reason the body_key characters are all \x00
                    # so we will use idx_body for the dict key
                            
                    vals, keys = next(lines).split("|")
                    vals = parse_floats(vals)
                    keys = keys.split(':')[-1].strip()
                    keys = keys.split(' ')
                    
                    current_case["outputs"]["body_forces"][idx_body] = {}
                
                    for k, v in zip(keys, vals):
                        current_case["outputs"]["body_forces"][idx_body][k] = v

            elif "s>  SURF" in line:
                # skip version line
                # and reference info
                for _ in range(4):
                    next(lines)                # skip a
                
                nsurfs, _ = next(lines).split("|")
                nsurfs = int(nsurfs)
                
                for idx_surf in range(1,nsurfs+1):
                    # skip next  line
                    for _ in range(1):
                        next(lines)
                    
                    surf_key = next(lines).strip()
                    
                    current_case["outputs"]["surface_forces"][surf_key] = {}
                    
                    # two data lines 
                    for _ in range(2):
                        vals, keys = next(lines).split("|")
                        vals = parse_floats(vals)
                        keys = keys.split(':')[-1].strip()
                        keys = keys.split(' ')
                        
                        for k, v in zip(keys, vals):
                            current_case["outputs"]["surface_forces"][surf_key][k] = v

            elif "s>  STRP" in line:
                # skip version line
                # and reference info
                for _ in range(5):
                    next(lines)                # skip a
                
                nsurfs, _ = next(lines).split("|")
                nsurfs = int(nsurfs)
                
                for idx_surf in range(1,nsurfs+1):
                    # skip next  line
                    for _ in range(1):
                        next(lines)
                    
                    surf_key = next(lines).strip()
                    
                    current_case["outputs"]["strip_forces"][surf_key] = {}
                    
                    surf_parms = next(lines).split("|")[0]
                    surf_parms = surf_parms.strip().split("   ")
                    # idx_surf = surf_parms[0]
                    nchord = int(surf_parms[1])
                    nspan = int(surf_parms[2])
                    first_strip = int(surf_parms[3])
                    
                    
                    # skip the surface data
                    for _ in range(4):
                        next(lines)
                    
                    # get the keys of the strips 
                    strip_vars = next(lines).strip().split(", ")  # skip header    
                    
                    # construct arrays of the strip data with the size of the nspan strips
                    for var in strip_vars:
                        var = var.strip()
                        current_case["outputs"]["strip_forces"][surf_key][var] = [0.0]*nspan
                    
                    # looop over the rows of data for each strip
                    # this is where you stopped
                    for idx_strip in range(nspan):
                        vals = parse_floats(next(lines))
                        
                        for idx, var in enumerate(strip_vars):
                            var = var.strip()
                            current_case["outputs"]["strip_forces"][surf_key][var][idx_strip] = vals[idx]
                    

            # --- Stability derivatives ---
            elif "Stability-axis derivatives" in line or "Geometry-axis derivatives" in line:
                
                if "Stability-axis derivatives" in line:
                    dict_key = "stability_derivatives"
                elif "Geometry-axis derivatives" in line:
                    dict_key = "body_axis_derivatives"
                
                # parse fixed blocks
                for _ in range(2):
                    next(lines)
                    for _ in range(6):
                        vals, keys = next(lines).split("|")
                        vals = parse_floats(vals)
                        keys = keys.split(':')[-1]
                        keys = keys.split(',')
                        
                        for k, v in zip(keys, vals):
                            current_case["outputs"][dict_key][k.strip()] = v
                
                # now do he control variables
                ncontrols, _ = next(lines).split("|")
                ncontrols = int(ncontrols)
                control_names = []
                for _ in range(ncontrols):
                    control_names.append(next(lines).strip())
                
                # read the matrix
                niters = 8 if dict_key == "stability_derivatives" else 6
                for _ in range(niters):
                    vals, keys = next(lines).split("|")
                    vals = parse_floats(vals)
                    key = keys.split(':')[-1]
                    
                    for idx, v in enumerate(vals):
                        k = key[:-1] + control_names[idx]
                        current_case["outputs"][dict_key][k.strip()] = v
                
                ndesign, _ = next(lines).split("|")
                ndesign = int(ndesign)
                if ndesign:
                    for _ in range(ndesign):
                        # skip dv name 
                        next(lines)
                        
                        # skip all the dv force and mom derivs
                        for _ in range(6):
                            next(lines)
                        
                        # skip CDffg and eg derivs
                        for _ in range(2):
                            next(lines)
                
                if dict_key == "stability_derivatives":
                    
                    # get the last two functions
                    val, _ = next(lines).split("|")
                    current_case["outputs"][dict_key]["Xnp"] = float(val)
                    val, _ = next(lines).split("|")
                    current_case["outputs"][dict_key]["spiral parameter"] = float(val)
                    
                
        # --- Body-axis derivatives ---

        # append last case
        if current_case:
            cases.append(current_case)
        

    return cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse AVL output into JSON.")
    parser.add_argument("--input_file", default="avl_output.txt", help="Path to AVL .txt output file")
    parser.add_argument("--output_file",default="test.json", help="Path to save parsed JSON")
    args = parser.parse_args()

    # with open(args.input_file, "r") as f:
        # file_text = f.read()

    cases = parse_avl_output(args.input_file)

    with open(args.output_file, "w") as f:
        json.dump(cases, f, indent=2)
