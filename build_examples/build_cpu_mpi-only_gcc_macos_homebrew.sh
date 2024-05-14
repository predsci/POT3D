#!/bin/bash
###################################################################
# This build assumes that you have an "mpif90" in your PATH 
#  set up to use your chosen MPI library and compiler.
##################################################################
#################################################################
# Please set the location of HDF5 include/library files and
#  the linker flags to match your installed version.
#
# Note! The HDF5 library needs to have been compiled with
#  the same compiler being used here and is loaded in the run-time
#  environment (e.g. LD_LIBRARY_PATH).
#################################################################

# Location of local hdf5 installed with same compiler being used for POT3D:
HDF5_INCLUDE_DIR="/usr/local/include/"
HDF5_LIB_DIR="/usr/local/lib/"
# Fortran HDF5 library flags (these can be version dependent):
HDF5_LIB_FLAGS="-lhdf5_fortran -lhdf5_hl_fortran -lhdf5 -lhdf5_hl"

###########################################################################
# Please set the compile flags based on your compiler and hardware setup.
###########################################################################

FFLAGS="-O3 -march=native"

###########################################################################
# If using NV HPC SDK for GPUs, with CUDA version >= 11.3, you can set 
# POT3D_CUSPARSE to "1" to link the cuSparse library, allowing you to set
# 'ifprec=2' in 'pot3d.dat' to yield ~2x speed improvement! 
# Warning!  Using ifprec=2 takes much more GPU memory than ifprec=1.
# You must also set CCFLAGS to use OpenACC in the C code.
###########################################################################

POT3D_CUSPARSE=0
CCFLAGS="-O3"

###########################################################################
###########################################################################
###########################################################################

POT3D_HOME=$PWD

cd ${POT3D_HOME}/src
echo "Making copy of Makefile..."
cp Makefile.template Makefile
echo "Modifying Makefile to chosen flags..."
sed -i "s#<FFLAGS>#${FFLAGS}#g" Makefile
sed -i "s#<CCFLAGS>#${CCFLAGS}#g" Makefile
sed -i "s#<POT3D_CUSPARSE>#${POT3D_CUSPARSE}#g" Makefile
sed -i "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" Makefile
echo "Building POT3D...."
make 1>build.log 2>build.err
echo "Copying POT3D executable from SRC to BIN..."
cp ${POT3D_HOME}/src/pot3d ${POT3D_HOME}/bin/pot3d
echo "Done!"

