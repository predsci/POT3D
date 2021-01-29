#!/bin/bash

###################################################################
# This build assumes that you have an "mpif90" in your PATH that is
# set up to use your chosen MPI library and compiler.
##################################################################

#################################################################
# Please set the location of your HDF5 include and library files. 
# Make sure the HDF5 library is compiled with 
# the same compiler currently being used and that the 
# library is in your run-time environment (e.g. LD_LIBRARY_PATH).
#################################################################

#HDF5_INCLUDE_DIR="/usr/include/hdf5/serial"
#HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"

HDF5_INCLUDE_DIR="/opt/psi/gnu/ext_deps/deps/hdf5/include"
HDF5_LIB_DIR="/opt/psi/gnu/ext_deps/deps/hdf5/lib"

###########################################################################
# Please set the HDF5 linker flags to match your installed version of hdf5.
###########################################################################

#HDF5_LIB_FLAGS="-lhdf5_serial_fortran -lhdf5_serialhl_fortran -lhdf5_serial -lhdf5_serial_hl"
HDF5_LIB_FLAGS="-lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl"

###########################################################################
# Please set the compile flags based on your compiler and hardware setup.
# Examples:
#   GNU (CPU):     FFLAGS="-O3 -mtune=native "
#   NV/PGI (CPU):  FFLAGS="-O3"
#   NV/PGI (GPU):  FFLAGS="-O3 -acc=gpu -gpu=cc60,cc70,cuda11.1 -Minfo=accel"
#   INTEL (CPU):   FFLAGS="-O3 -fp-model precise -assume byterecl -heap-arrays -xCORE-AVX2 -axCORE-AVX512"
###########################################################################

FFLAGS="-O3 -mtune=native"

###########################################################################
###########################################################################
###########################################################################

POT3D_HOME=$PWD

cd ${POT3D_HOME}/src
cp Makefile.template Makefile
sed -i "s#<FFLAGS>#${FFLAGS}#g" Makefile
sed -i "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" Makefile
make

echo "cp ${POT3D_HOME}/src/pot3d ${POT3D_HOME}/bin/pot3d"

cp ${POT3D_HOME}/src/pot3d ${POT3D_HOME}/bin/pot3d

