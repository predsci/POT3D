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

#  Ubuntu 20.x:
HDF5_INCLUDE_DIR="/usr/include/hdf5/serial"
HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"
HDF5_LIB_FLAGS="-lhdf5_serial_fortran -lhdf5_serialhl_fortran -lhdf5_serial -lhdf5_serial_hl"

#  Locally installed older version example:
#HDF5_INCLUDE_DIR="/opt/psi/nv/ext_deps/deps/hdf5/include"
#HDF5_LIB_DIR="/opt/psi/nv/ext_deps/deps/hdf5/lib"
#HDF5_LIB_FLAGS="-lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl"

###########################################################################
# Please set the compile flags based on your compiler and hardware setup.
###########################################################################

FFLAGS="-O3"

###########################################################################
# Examples:
#   GCC              (CPU MPI only):    FFLAGS="-O3 -march=native"
#   GCC              (CPU MPI+threads): FFLAGS="-O3 -march=native 
#                                               -ftree-parallelize-loops=${OMP_NUM_THREADS}"
#   NVIDIA HPC SDK   (CPU MPI only):    FFLAGS="-O3 -march=native"
#   NVIDIA HPC SDK   (CPU MPI+threads): FFLAGS="-O3 -march=native 
#                                               -stdpar=multicore -acc=multicore"
#   NVIDIA HPC SDK   (GPU MPI+GPU):     FFLAGS="-O3 -march=native 
#                                               -stdpar=gpu -acc=gpu -Minfo=accel 
#                                               -gpu=cc80,cuda11.6,nomanaged"
#   INTEL HPC SDK    (CPU MPI only):    FFLAGS="-O3 -xHost -assume byterecl 
#                                               -heap-arrays"
#   INTEL HPC SDK    (CPU MPI+threads): FFLAGS="-O3 -xHost -assume byterecl 
#                                               -heap-arrays -mp"
###########################################################################

###########################################################################
# If using NV HPC SDK for GPUs, with CUDA version >= 11.3, you can set 
#  the following to "1" to link the cuSparse library, allowing you to set
#  'ifprec=2' in 'pot3d.dat' to yield ~2x speed improvement! 
# Warning!  Using ifprec=2 takes much more GPU memory than ifprec=1.
###########################################################################

POT3D_CUSPARSE=0

###########################################################################
###########################################################################
###########################################################################

POT3D_HOME=$PWD

cd ${POT3D_HOME}/src
cp Makefile.template Makefile
sed -i "s#<FFLAGS>#${FFLAGS}#g" Makefile
sed -i "s#<POT3D_CUSPARSE>#${POT3D_CUSPARSE}#g" Makefile
sed -i "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" Makefile
echo "make 1>build.log 2>build.err"
make 1>build.log 2>build.err

echo "cp ${POT3D_HOME}/src/pot3d ${POT3D_HOME}/bin/pot3d"

cp ${POT3D_HOME}/src/pot3d ${POT3D_HOME}/bin/pot3d

