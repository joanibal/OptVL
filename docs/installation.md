# Installation Guide for OptVL

## Installing with pip
The recommended method to install `optvl` is via pip:

```shell
pip install optvl
```
This package is packaged with an OpenBLAS linear solver for quicker analysis.

### Supported Platforms
Currently, `optvl` supports Linux, macOS (Apple Silicon and Intel), and Windows!

## Building Locally
If you would like to build `optvl` manually, follow the steps below:

1. Clone the repository to your local machine.
2. Navigate to the root directory and run:
   ```
   pip install .
   ```
You will need the following dependencies
1. meson
2. meson-python
3. ninja

## Post-processing
You can use `matplotlib` for viewing the geometry and 3D surface CP data. 
You can also install it with `pip`. 
You can also write out 3D data in the ASCII tecplot format that can be opened in either Tecplot or Paraview. 

See the section on post-processing for more details about how each of these are used. 

