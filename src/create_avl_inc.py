#!/usr/bin/env python


import sys


def replace_text(input_file, output_file, search_text, replace_text):
    with open(input_file, "r") as file:
        file_contents = file.read()

    # Replace all occurrences of the search_text with replace_text
    updated_contents = file_contents.replace(search_text, replace_text)

    with open(output_file, "w") as file:
        file.write(updated_contents)


if __name__ == "__main__":
    # Arguments: input_file output_file search_text replace_text

    import os

    base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder

    input_file = os.path.join(base_dir, "includes/AVL.INC.in")
    output_file = os.path.join(base_dir, "includes/AVL.INC")
    search_text = "@OSNVMAX@"

    import platform

    if platform.system() == "Darwin":
        nvmax = 6000
    elif platform.system() == "Linux":
        nvmax = 6000
    elif platform.system() == "Windows":
        nvmax = 5000
    else:
        raise NotImplementedError("system platform not found")

    replacement_text = str(nvmax)

    replace_text(input_file, output_file, search_text, replacement_text)
