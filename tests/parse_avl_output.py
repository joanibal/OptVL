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
                for _ in range(2):
                    vals, keys = next(lines).split("|")
                    vals = parse_floats(vals)
                    keys = keys.split(',')
                    
                    for k, v in zip(keys, vals):
                        current_case["outputs"]["reference_data"][k.strip()] = v

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
