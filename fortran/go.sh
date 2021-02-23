#!/bin/bash
# rm -f harms *.mod
# #DBGFLAG="-ffree-line-length-none -g -fbacktrace -fcheck=all -pedantic -Wall -Wextra -W -Wno-unused-function -fopenmp"
# # OPTFLAG="-ffree-line-length-none -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe"
# DBGFLAG="-O0 -g -traceback -check all -check bounds -debug all -gen-interfaces -warn interfaces"
# OPTFLAG="-O3 -funroll-loops -fomit-frame-pointer "
# PRFFLAG="-profile-functions -profile-loops=all -profile-loops-report=2"
# FLAG=$OPTFLAG
# echo ifort $FLAG  -o harms UTIL.f90 HARMONICS.f90 SHAPE_noali.f90 main.f90 FFts/fftpack5.o -L./lib/lapack-3.9.0/ -llapack -lrefblas
# ifort $FLAG  -o harms UTIL.f90 HARMONICS.f90 SHAPE_noali.f90 main.f90 FFts/fftpack5.o -L./lib/lapack-3.9.0/ -llapack -lrefblas
# echo "Compiled"
# rm -f *.o
rm -f harms
make dbg=0 prof=0
./harms noBgrint.mfs