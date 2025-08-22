# Config File for LINUX and GFORTRAN Compiler
AR       = ar
AR_FLAGS = -rvs

FF90       = gfortran



FF90_FLAGS = -fdefault-real-8 -fdefault-double-8 -O2 -fPIC -Wno-align-commons -std=legacy -C -mmacosx-version-min=13.6 
# FF90_FLAGS = -fdefault-real-8 -fdefault -double-8 -O0 -fPIC -Wno-align-commons -Werror=line-truncation -std=legacy -C -mmacosx-version-min=13.6 -g -fcheck=bounds -finit-real=snan -finit-integer=-999999 -ftrapping-math   -ftrapv


C_FLAGS = -O2 -fPIC -mmacosx-version-min=13.6

F2PY = f2py
F2PY_FF90 = gfortran

PYTHON = python

LINKER_FLAGS =  -undefined dynamic_lookup
  