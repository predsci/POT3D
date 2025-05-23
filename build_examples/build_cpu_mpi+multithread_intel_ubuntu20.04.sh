#!/bin/bash
#################################################################
# Please set the location of HDF5 include/library files and
#  the linker flags to match your installed version.
#
# Build script template for POT3D.
#
# See the build scripts in the "/build_examples" 
# folder for examples of what to insert here.
# Please do not edit this file, rather copy it to a new 
# file with your own configuration options inserted.
# If your configuration is unique, feel free to add it to the 
# /build_examples folder.
# Note that your build script must be launched from 
# this lowest level directory (not from the build_examples folder).
#
#################################################################
#################################################################
# Enter your MPI compiler (typically "mpif90").
#################################################################

FC="mpif90 -f90=ifx"

#################################################################
# Please set the location of the HDF5 include & library files. 
# Make sure the HDF5 LIBRARY is COMPILED with 
# the SAME COMPILER used here, and is in the run-time environment.
#################################################################

HDF5_INCLUDE_DIR="/usr/include/hdf5/serial"
HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"

##################################################################
# Please set the HDF5 linker flags to match the installed version.
##################################################################

HDF5_LIB_FLAGS="-lhdf5_serial_fortran -lhdf5_serialhl_fortran -lhdf5_serial -lhdf5_serial_hl"

###########################################################################
# Please set the compile flags based on your compiler and hardware setup.
###########################################################################

FFLAGS="-O3 -xHost -fp-model precise -heap-arrays -fopenmp-target-do-concurrent -fiopenmp -fopenmp-targets=spir64_x86_64"

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

cX="\033[0m"
cR="\033[1;31m"
cB="\033[1;34m"
cG="\033[32m"
cC="\033[1;96m"
cM="\033[35m"
cY="\033[1;93m"
Bl="\033[1;5;96m"
echo="echo -e"

${echo} "${cG}=== STARTING POT3D BUILD ===${cX}"
${echo} "==> Entering src directory..."
pushd ${POT3D_HOME}/src > /dev/null
${echo} "==> Removing old Makefile..."
if [ -e Makefile ]; then
  \rm Makefile
fi 
${echo} "==> Generating Makefile from Makefile.template..."
sed \
  -e "s#<FC>#${FC}#g" \
  -e "s#<FFLAGS>#${FFLAGS}#g" \
  -e "s#<CCFLAGS>#${CCFLAGS}#g" \
  -e "s#<POT3D_CUSPARSE>#${POT3D_CUSPARSE}#g" \
  -e "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" \
  -e "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" \
  -e "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" \
  Makefile.template > Makefile
${echo} "==> Compiling code..."
make clean 1>/dev/null 2>/dev/null ; make 1>build.log 2>build.err
if [ ! -e pot3d ]; then
  ${echo} "${cR}!!> ERROR!  pot3d executable not found.  Build most likely failed."
  ${echo} "            Contents of src/build.err:"
  cat build.err
  ${echo} "${cX}"
  exit 1
fi
${echo} "==> Copying pot3d executable to: ${POT3D_HOME}/bin/pot3d"
cp pot3d ${POT3D_HOME}/bin/pot3d
${echo} "${cG}==> Build complete!${cX}"
${echo}      "    Please add the following to your shell startup (e.g. .bashrc, .profile, etc.):"
${echo} "${cC}    export PATH=${POT3D_HOME}/bin:\$PATH${cX}"

