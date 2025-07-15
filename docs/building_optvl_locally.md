# Building OptVL Locally
If you are developing OptVL you will want to install OptVL locally.
After cloning the repo and navigating to the root directory there are two methods to actually build the python package, the meson way or the make way. 

## Using Meson
This method exists to integrate with the tools necessary to package the binaries into wheels for distribution on PYPI, 
The steps are 
1. install meson and ninja
2. run `pip install -e .`

The `pyproject.toml` should hold all the information pip needs to use the meson backend for compiling the code. 
To modify aspects of the build, such as the compilation flags, see the `meson.build` file. 

## Using Make
This method is more convenient for quick compilation during editing. The steps are

1. Copy a config file from `config/defaults` to the config directory, i.e. `cp  config/defaults/config.LINUX_GFORTRAN.mk config/config.mk`
2. Create the shared library file with `make`
3. Install the python wrapper with `python setup_deprecated.py develop`

After modifying the Fortran code only `make` needs to be run again. 
If you only modify the python code then there is no need to rerun anything. 

