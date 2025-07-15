# Common Development Task
Below is some information on a few development tasks that you may be interested in. 

## Modify the python wrapper or OpenMDAO wrapper.
The file `optvl_class.py` wraps the Fortran object file and is where the user API is defined. 
For most development tasks this is what you will need to modify. 
The file `om_wrapper.py` wraps the `OVLSolver` class for use with OpenMDAO. 
If you want to change something about the OpenMDAO interface, you'll need to start here. 

## Run AD on modified code
Tapenade 3.16 is used to generate the derivative routines. 
After installing Tapenade (only really works on Linux), navigate to the `/src/ad_src`.
From here run `make -f Makefile_tapenade`.
This will make both the forward and reverse routines by default. 
After Tapenade generates the routines a few loop dimensions will need to be changed for speed. 
I can do this step for you if you make a pull request.
But in general there are a few places where we can loop over NVOR instead of NVMAX.


## Expose an additional variable from the Fortran level
1. in `AVL.INC.in` add it to the COMMON line AND declare its type above
2. replace instances of the local variable in the code with the global variable from the COMMON block
3. navigate to `src/includes` and run `gen_ad_inc.py` 
4. rebuild the library 

## Releasing a new version
1. Bump the version number in `pyproject.toml` AND in `meson.build`.
    - I cannot figure out a way to single source this so it will have to be done in two steps for now.
2. Then create a new release on GitHub. Create a new tagged version as part of the release
3. Check that that new version has been released on PyPI as intended
