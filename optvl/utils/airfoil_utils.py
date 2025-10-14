import numpy as np
import re

def read_coordinates_files(filename, headerlines=0):
    """
    Reads a '.dat' style airfoil coordinate file,
    with each coordinate on a new line and each
    line containing an xy pair separate by whitespace.

    Args:
        filename : str
            dat file from which to read data
        headerlines : int
            the number of lines to skip at the beginning of the file to reach the coordinates

    Returns:
        X : Ndarray [N,2]
            The coordinates read from the file
    """
    with open(filename, "r") as f:
        for _i in range(headerlines):
            f.readline()
        r = []
        while True:
            line = f.readline()
            if not line:
                break  # end of file
            if line.isspace():
                break  # blank line
            if re.search("[a-zA-Z]", line):  # line contains chars
                continue
            r.append([float(s) for s in line.split()])

            X = np.array(r)

    return X
