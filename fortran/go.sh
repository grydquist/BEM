#!/bin/bash
rm -f harms *.mod
#DBGFLAG="-ffree-line-length-none -g -fbacktrace -fcheck=all -pedantic -Wall -Wextra -W -Wno-unused-function -fopenmp"
# OPTFLAG="-ffree-line-length-none -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe"
DBGFLAG="-O0 -g -traceback -check all -check bounds -debug all -gen-interfaces -warn interfaces"
OPTFLAG="-O3 -funroll-loops -fomit-frame-pointer "
# set -x
ifort $DBGFLAG  -o harms UTIL.f90 Fast_Fourier.f90 HARMONICS.f90 SHAPE.f90 main.f90 -L./lib/lapack-3.9.0/ -llapack -lrefblas
echo "Compiled"
rm -f *.o
./harms grint.mfs