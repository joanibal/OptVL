import os
import json
import subprocess
from generate_avl_scripts import generate_avl_script
from parse_avl_output import parse_avl_output

abs_path = os.path.abspath(__file__)
test_dir = os.path.dirname(abs_path)
avl_exe = os.path.join(test_dir, "../external/avl_db")

# this data is used to generate the avl scripts
# fmt:off
case_data = {
    "unconstrained_supra": {
        "geometry_file": "../geom_files/supra.avl",
        "mass_file": "../geom_files/supra.mass",
        "alphas":      [0, 45, 0, -45],
        "betas":       [0, 0,  0, -45],
        "roll_rates":  [0, 0,  0, 0.1],
        "pitch_rates": [0, 0,  0, 0.2],
        "yaw_rates":   [0, 0,  0, 0.3],
        "controls" : {
            "D1":      [0, 0,  10, 10],
            "D2":      [0, 0,  10, 10],
            "D3":      [0, 0,  10, 10],
        },
    },
    "constrained_supra": {
        "geometry_file": "../geom_files/supra.avl",
        "mass_file": "../geom_files/supra.mass",
        "CLs": [1.0, 1.1, 1.2],
        "betas": [0, 10, -12],
        "roll_rates": [0, 0, 0.1],
        "pitch_rates": [0, 0, 0.2],
        "yaw_rates": [0, 0, 0.3],
        "pre_commands": [("D2", "RM", 0), ("D3", "PM", 0), ("D4", "YM", 0)],
    },
    "unconstrained_aircraft": {
        "geometry_file": "../geom_files/aircraft.avl",
        "mass_file": "../geom_files/aircraft.mass",
        "alphas":      [0, 45,  0,   0,   0,   0,  0,  0,  -45],
        "betas":       [0,  0, 45,   0,   0,   0,  0,  0,  -45],
        "roll_rates":  [0,  0,  0, 0.1,   0,   0,  0,  0,  0.1],
        "pitch_rates": [0,  0,  0,   0, 0.1,   0,  0,  0, -0.1],
        "yaw_rates":   [0,  0,  0,   0,   0, 0.1,  0,  0,  0.1],
        "controls" : {
            "D1":      [0,  0,  0,   0,   0,   0, 10,  0,  10],
            "D2":      [0,  0,  0,   0,   0,   0,  0, 10,  10],
        },
    },
    "constrained_aircraft": {
        "geometry_file": "../geom_files/aircraft.avl",
        "mass_file": "../geom_files/aircraft.mass",
        "CLs": [0.0, 0.5, 1.01],
        "betas": [0, 10, -12],
        "roll_rates":  [0, 0, 0.5],
        "pitch_rates": [0, 0, 0.33],
        "yaw_rates":   [0, 0, 0.01],
        "pre_commands": [("D1", "PM", 0), ("D2", "YM", 0)],
    },
}
# fmt:on
for case in case_data:
    avl_script_file = f"avl_analysis_scripts/{case}.txt"
    avl_script_file = os.path.join(test_dir, f"avl_analysis_scripts/{case}.txt")
    # generate the script first
    _ = generate_avl_script(**case_data[case], output_file=avl_script_file)

    # run the script
    avl_run_file = os.path.join(test_dir, f"avl_analysis_references/{case}.output")
    cmd = [avl_exe, "<", avl_script_file, "| tee", avl_run_file]

    with open(avl_script_file, "r") as fin, open(avl_run_file, "w") as fout:
        result = subprocess.run(
            [avl_exe], stdin=fin, stdout=fout, stderr=subprocess.PIPE, cwd=test_dir, text=True, timeout=120
        )
    print(f"{case} piped through avl and exited with:", result.returncode)

    cases = parse_avl_output(avl_run_file)
    reference_file = os.path.join(test_dir, f"avl_analysis_references/{case}.json")
    with open(reference_file, "w") as fout:
        json.dump(cases, fout, indent=2)
