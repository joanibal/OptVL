"""
Common utilities for writing Python dictionaries and arrays to .py files
"""

import numpy as np


def write_value(value, indent_level=0):
    """
    Convert a value to its string representation for writing to file.

    Handles numpy arrays, numpy scalars, dicts, lists, and primitives.

    Parameters
    ----------
    value : any
        The value to convert to a string representation
    indent_level : int
        The indentation level (each level = 4 spaces)

    Returns
    -------
    str
        String representation suitable for writing to a .py file
    """
    indent = "    " * indent_level

    if isinstance(value, np.ndarray):
        return f"np.array({value.tolist()}, dtype=np.{value.dtype})"
    elif isinstance(value, (np.int32, np.int64)):
        return f"np.int32({int(value)})"
    elif isinstance(value, (np.float32, np.float64)):
        return f"np.float64({float(value)})"
    elif isinstance(value, str):
        return f'"{value}"'
    elif isinstance(value, dict):
        lines = ["{\n"]
        for k, v in value.items():
            lines.append(f'{indent}    "{k}": {write_value(v, indent_level + 1)},\n')
        lines.append(f'{indent}}}')
        return ''.join(lines)
    elif isinstance(value, (list, tuple)):
        return repr(value)
    elif isinstance(value, (int, float)):
        return repr(value)
    else:
        return repr(value)


def write_dict_to_py(output_file, data_dict, var_name='input_dict', description=''):
    """
    Write a dictionary to a .py file.

    Parameters
    ----------
    output_file : str
        Path to the output .py file
    data_dict : dict
        Dictionary to write
    var_name : str
        Variable name to use in the output file (default: 'input_dict')
    description : str
        Description to include in the docstring
    """
    with open(output_file, 'w') as f:
        # Write docstring
        f.write(f'"""Auto-generated geometry file\n')
        if description:
            f.write(f'{description}\n')
        f.write('"""\n')
        f.write('import numpy as np\n\n')

        # Write the dictionary
        f.write(f'{var_name} = {{\n')
        for key, value in data_dict.items():
            f.write(f'    "{key}": {write_value(value, 1)},\n')
        f.write('}\n')


def write_array_to_py(output_file, array, var_name='mesh', description=''):
    """
    Write a numpy array to a .py file.

    Parameters
    ----------
    output_file : str
        Path to the output .py file
    array : np.ndarray
        Array to write
    var_name : str
        Variable name to use in the output file (default: 'mesh')
    description : str
        Description to include in the docstring
    """
    with open(output_file, 'w') as f:
        # Write docstring
        f.write(f'"""Auto-generated array file\n')
        if description:
            f.write(f'{description}\n')
        f.write('"""\n')
        f.write('import numpy as np\n\n')

        # Write the array
        f.write(f'{var_name} = np.array({array.tolist()}, dtype=np.{array.dtype})\n')
