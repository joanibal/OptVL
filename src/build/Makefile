# Include user supplied configuration file
include ../../config/config.mk

# Group Fortran compiler flags
FF90_ALL_FLAGS = $(FF90_FLAGS)

CC_ALL_FLAGS     = $(C_FLAGS)
# Include full list of files
include fileList

# Include full list of directories
include directoryList

# Include the list of rules
include rules

# Set the make VPATH variable to the "dirs" variable from
# directorylist. We must first append the '../'
dirs:=$(addprefix ../,$(dirs))
VPATH:=$(dirs)

# Add ../ prefix to all source files, remove directories
f90Files:=$(addprefix ../,$(f90Files))
cFiles:=$(addprefix ../,$(cFiles))
f77Files:=$(addprefix ../,$(f77Files))

f90FilesNoDir=$(notdir $(f90Files))
cFilesNoDir=$(notdir $(cFiles))
f77FilesNoDir=$(notdir $(f77Files))

# Generate two separate list of .F90 and .f90 files using the filter command

# Finally convert all source files to .o
OFILES=$(f90FilesNoDir:%.f90=%.o) $(f77FilesNoDir:%.f=%.o) $(cFilesNoDir:%.c=%.o)

# Compile sources


default: lib ../f2py/libavl.pyf
# Generate Python inlude directory
	 $(eval PYTHON_INCLUDES = $(shell python-config --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print(numpy.get_include())'))
	@echo "#------------------------------------------------------#"
	@echo Numpy Include Directory: $(NUMPY_INCLUDES)
	@echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) ../f2py/get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"


# take all the of f2py steps one at a time to avoid using a backend 
# Run f2py to get libavl-f2pywrappers.f and libavlmodule.c
	$(F2PY) ../f2py/libavl.pyf
	
# Compile c wrapper. Don't use CC_ALL_FLAGS...PETSc wil F-up this command.
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libavlmodule.c

# Compile fortranobject needed by all f2py modules
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compiled f2py-generated wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I../includes -I./ -c libavl-f2pywrappers.f

# Final Link:
	$(FF90) $(LINKER_FLAGS) -shared fortranobject.o libavlmodule.o libavl-f2pywrappers.o libavl.a -o libavl.so
	$(PYTHON) importTest.py libavl
	mv libavl.*so ../../optvl/libavl.so

# # asetup uses avl_heap_inc, so make sure .mod exists
# asetup.o: avl_heap_inc.mod

# build the module first
avl_heap_inc.o avl_heap_inc.mod: ../avl_heap_inc.f90
	$(FF90) $(FF90_ALL_FLAGS) -cpp -I ../includes -c $<  -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

avl_heap_diff_inc.o avl_heap_diff_inc.mod: ../avl_heap_diff_inc.f90
	$(FF90) $(FF90_ALL_FLAGS) -cpp -I ../includes -c  $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

sources: $(OFILES) 

lib: sources
	$(AR) $(AR_FLAGS) libavl.a $(OFILES) 
