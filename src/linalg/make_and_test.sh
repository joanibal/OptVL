#!/bin/sh

# Set compiler flags
FFLAGS="-O0 -g -pedantic-errors -fcheck=bounds"
# echo "Compiling with flags: $FFLAGS"

# Compile source files
gfortran $FFLAGS -c linalg.f90
gfortran $FFLAGS -c test_linalg.f90

# Link object files and create executable
gfortran $FFLAGS -o test_linalg linalg.o test_linalg.o -llapack -lblas

# Run the executable
./test_linalg