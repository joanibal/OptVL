def generate_avl_script(
    geometry_file,
    mass_file=None,
    alphas=None,
    CLs=None,
    betas=None,
    roll_rates=None,
    pitch_rates=None,
    yaw_rates=None,
    controls=None,
    pre_commands=None,
    output_file="avl_script.txt",
):
    """
    Generate an AVL script and save it to a file.

    Parameters:
    ------------
    geometry_file : str
        Path to the AVL geometry file.
    alphas : list[float]
        List of angles of attack to sweep through.
    betas : list[float]
        List of betas to sweep through.
    commands : list[tuple], optional
        List of tuples (e.g. [("D1", "PM", 0), ("D2", "RM", 0)])
        that are inserted before the angle sweeps.
    mass_file : str, optional
        Path to the AVL mass file. If None, skip mass block.
    output_file : str, optional
        Path to save the generated script (default: "avl_script.txt").

    Returns:
    --------
    str : The generated AVL script as a string.
    """

    lines = []

    # Geometry file
    lines.append(f"load {geometry_file}")

    # Optional mass block
    if mass_file:
        lines.append(f"mass {mass_file}")
        lines.append("mset 0")

    # Enter OPER menu
    lines.append("oper")
    lines.append("MRF")

    # Insert user-defined commands
    if pre_commands:
        for cmd in pre_commands:
            # Each command is a tuple like ("D1", "PM", 0)
            lines.append(" ".join(map(str, cmd)))

    if CLs is None:
        # Angle of attack sweeps
        for i in range(len(alphas)):
            lines.append(f"a a {alphas[i]}")
            lines.append(f"b b {betas[i]}")
            lines.append(f"r r {roll_rates[i]}")
            lines.append(f"p p {pitch_rates[i]}")
            lines.append(f"y y {yaw_rates[i]}")
            if controls is not None:
                for con in controls:
                    lines.append(f"{con} {con} {controls[con][i]}")
            lines.append("x")
            lines.append("FT")
            lines.append("")
            lines.append("FB")
            lines.append("")
            lines.append("FN")
            lines.append("")
            lines.append("ST")
            lines.append("")
            lines.append("SB")
            lines.append("")
            lines.append("FS")
            lines.append("")

    if alphas is None:
        # CL sweeps
        for i in range(len(CLs)):
            lines.append(f"a c {CLs[i]}")
            lines.append(f"b b {betas[i]}")
            lines.append(f"r r {roll_rates[i]}")
            lines.append(f"p p {pitch_rates[i]}")
            lines.append(f"y y {yaw_rates[i]}")
            if controls is not None:
                for con in controls:
                    lines.append(f"{con} {con} {controls[con][i]}")
            lines.append("x")
            lines.append("FT")
            lines.append("")
            lines.append("FB")
            lines.append("")
            lines.append("FN")
            lines.append("")
            lines.append("ST")
            lines.append("")
            lines.append("SB")
            lines.append("")
            lines.append("FS")
            lines.append("")

    # Quit AVL
    lines.append("")
    lines.append("")
    lines.append("")
    lines.append("q")
    lines.append("q")

    script_content = "\n".join(lines)

    # Write to file
    with open(output_file, "w") as f:
        f.write(script_content)

    return script_content


if __name__ == "__main__":
    # Example usage:
    script = generate_avl_script(
        geometry_file="./geom_files/supra.avl",
        mass_file="./geom_files/supra.mass",
        alphas=[0, 45, -45],
        betas=[0, 0, -45],
        output_file="avl_analysis_scripts/unconstrained_supra.txt",
    )
    script = generate_avl_script(
        geometry_file="./geom_files/supra.avl",
        mass_file="./geom_files/supra.mass",
        # alphas=[0, 0, -10],
        CLs=[1.0, 1.1, 1.2],
        betas=[0, 10, -12],
        commands=[("D2", "RM", 0), ("D3", "PM", 0), ("D4", "YM", 0)],
        output_file="avl_analysis_scripts/constrained_supra.txt",
    )
