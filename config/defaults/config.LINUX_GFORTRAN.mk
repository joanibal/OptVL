# Config File for LINUX and GFORTRAN Compiler
AR       = ar
AR_FLAGS = -rvs

FF90       = gfortran
FF90_FLAGS = -O2 -fdefault-real-8 -fdefault-double-8 -fPIC -Wno-align-commons -std=legacy
# FF90_FLAGS = -O0 -fdefault-real-8 -fdefault-double-8 -fPIC -Wno-align-commons -std=legacy -g -fcheck=bounds -Werror=line-truncation

C_FLAGS = -O2 -fPIC -g
# C_FLAGS = -O0 -fPIC -g

F2PY = f2py
F2PY_FF90 = gfortran

PYTHON = python

LINKER_FLAGS =  