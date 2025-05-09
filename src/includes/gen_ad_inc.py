import re


def process_fortran_include_file(input_filename, output_filename, ad_ext="_DIFF"):
    variables = {}
    include_lines = []

    def get_line(idx_line):
        line = lines[idx_line]
        # remove fortan comments
        line = line.split("!")[0]
        # trim whitespace
        line = line.strip()
        return line

    ignore_blocks = ["UN_R", "TIME_R", "PLOT_R", "VRTX_S"]
    ignore_vars = ["VERSION", "DTR", "PI", "DTIMED", "AICN_LU(NVMAX,NVMAX)", "EXEC_TOL"]

    with open(input_filename, "r") as input_file:
        lines = input_file.readlines()

    idx_line = 0
    while idx_line < len(lines):
        line = get_line(idx_line)

        # search the line for keywords
        if "include" in line.lower():
            include_lines.append(line)

        elif "COMMON" in line:
            block_name = line.split("/")[1]

            if block_name in ignore_blocks:
                print("ignoring block", block_name)
                idx_line += 1
                continue

            # only add the variables if they are from a block of real variables
            # the file naming convention is that such block names end with an 'r'
            if block_name[-1].lower() == "r":
                variables[block_name] = []
                # loop until the end of the block
                idx_line += 1
                while True:
                    line = get_line(idx_line)
                    if line[0].lower() == "c":
                        idx_line += 1
                        continue
                    # remove & from the line
                    line = line.replace("& ", "")
                    # print(line)
                    is_last = line[-1] != ","

                    regex = r",(?![^(]*\))"  # Matches commas that are not within parentheses
                    _vars = re.split(regex, line)
                    
                    for _var in _vars:
                        # import pdb; pdb.set_trace()
                        _var = _var.strip()
                        if _var != "" and _var not in ignore_vars:
                            variables[block_name].append(_var)
                            
                        if _var in ignore_vars:
                            print(f"ignoring {_var}")
                        # import pdb; pdb.set_trace()

                    if is_last or idx_line >= len(lines) - 1:
                        break
                    else:
                        idx_line += 1

        idx_line += 1
    
    with open(output_filename, "w") as output_file:
        output_file.write("C    automatically generated by gen_ad_inc.py\n")
        output_file.write(f"C    must be precended by 'include AVL_ad_seed'\n")

        # add the include lines
        # for line in include_lines:
        #     output_file.write(f"      {line}\n")
        # output_file.write(f"\n\n")

        for block in variables:
            # declare the type of the variables in the block
            for idx_var, var in enumerate(variables[block]):
                if "(" in var:
                    # if size is give in variable, insert _diff inbetween variable and size
                    var, size = var.split("(")
                    new_var = f"{var}{ad_ext}"

                else:
                    new_var = var + ad_ext

                output_file.write(f"      real(kind=avl_real) {new_var}\n")

            output_file.write(f"      COMMON /{block}{ad_ext}/\n")
            for idx_var, var in enumerate(variables[block]):
                if idx_var == len(variables[block]) - 1:
                    is_last = True
                else:
                    is_last = False

                if "(" in var:
                    # if size is give in variable, insert _diff inbetween variable and size
                    var, size = var.split("(")
                    size = size.replace(")", "")
                    new_var = f"{var}{ad_ext}({size})"

                else:
                    new_var = var + ad_ext

                output_file.write(f"     & {new_var}{'' if is_last else ','}\n")

            output_file.write("C\n")


if __name__ == "__main__":
    process_fortran_include_file("AVL.INC.in", "AVL_ad_seeds.inc")
