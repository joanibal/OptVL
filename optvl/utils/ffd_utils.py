# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np

# =============================================================================
# Local modules
# =============================================================================


def write_FFD_file(mesh, mx, my, filename="ffd", cushion=0.05):
    """Utility functions that generates a box FFD around a wing mesh.
    Based on the utility function in OpenAeroStruct.

    Args:
        mesh: mesh array for the surface to be embedded in the FFD np.ndarray (nx,ny,3)
        mx: number of control points in the x-direction
        my: number of control points in the y-direction
        filename: name of the ffd file to write without extension.Defaults to "ffd.dat".
        cushion: Amount of cushion to apply from the LE/TE and tips. Defaults to 0.05.

    Returns:
        filename of the written FFD file as a str
    """
    nx, ny = mesh.shape[:2]

    half_ffd = np.zeros((mx, my, 3))

    LE = mesh[0, :, :]
    TE = mesh[-1, :, :]

    half_ffd[0, :, 0] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), LE[:, 0])
    half_ffd[0, :, 1] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), LE[:, 1])
    half_ffd[0, :, 2] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), LE[:, 2])

    half_ffd[-1, :, 0] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), TE[:, 0])
    half_ffd[-1, :, 1] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), TE[:, 1])
    half_ffd[-1, :, 2] = np.interp(np.linspace(0, 1, my), np.linspace(0, 1, ny), TE[:, 2])

    for i in range(my):
        half_ffd[:, i, 0] = np.linspace(half_ffd[0, i, 0], half_ffd[-1, i, 0], mx)
        half_ffd[:, i, 1] = np.linspace(half_ffd[0, i, 1], half_ffd[-1, i, 1], mx)
        half_ffd[:, i, 2] = np.linspace(half_ffd[0, i, 2], half_ffd[-1, i, 2], mx)

    half_ffd[0, :, 0] -= cushion
    half_ffd[-1, :, 0] += cushion
    half_ffd[:, 0, 1] -= cushion
    half_ffd[:, -1, 1] += cushion

    bottom_ffd = half_ffd.copy()
    bottom_ffd[:, :, 2] -= cushion

    top_ffd = half_ffd.copy()
    top_ffd[:, :, 2] += cushion

    ffd = np.vstack((bottom_ffd, top_ffd))

    filename = filename + ".xyz"

    with open(filename, "w") as f:
        f.write("1\n")
        f.write("{} {} {}\n".format(mx, 2, my))
        x = np.array_str(ffd[:, :, 0].flatten(order="F"))[1:-1] + "\n"
        y = np.array_str(ffd[:, :, 1].flatten(order="F"))[1:-1] + "\n"
        z = np.array_str(ffd[:, :, 2].flatten(order="F"))[1:-1] + "\n"

        f.write(x)
        f.write(y)
        f.write(z)

    return filename