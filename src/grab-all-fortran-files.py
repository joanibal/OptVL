#!/usr/bin/env python

# # Standard Python modules
import os


def parse_file(file_name):
    # Initialize lists for each file type
    f90Files = []
    f77Files = []
    cFiles = []

    # Define a dictionary to map identifiers to lists
    list_map = {
        "f90Files =": f90Files,
        "f77Files =": f77Files,
        "cFiles =": cFiles,
    }

    # Initialize current list as None
    current_list = None

    # Open and read the file line by line
    with open(file_name, "r") as file:
        for line in file:
            # Remove leading/trailing white space
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith("#") or not line:
                continue
            # Check if the line is a new list identifier
            new_list = False
            for key in list_map:
                if line.startswith(key):
                    current_list = list_map[key]
                    new_list = True
                    break

            # Otherwise, if a list is currently being filled, add the file to the list
            if current_list is not None and not new_list:
                # trim whitespace and trailing \
                line = line.strip()
                if line.endswith("\\"):
                    line = line[:-1]
                # prepend src/ to the file name
                line = "src/" + line

                # add to list
                current_list.append(line)

    # Print each file list
    # print("f90Files = ", f90Files)
    # print("f77Files = ", f77Files)
    # print("cFiles = ", cFiles)
    for file in f90Files:
        print(file)
    for file in f77Files:
        print(file)
    for file in cFiles:
        print(file)


# Call the function with the file name
parse_file("src/build/fileList")  # replace "filename.txt" with the name of your file
