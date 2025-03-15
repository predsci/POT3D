set -ex

if [[ "$(uname)" == "Linux" ]]; then
  CC=gcc
else
  CC=clang
fi

FC=lfortran
FFLAGS="--fast --skip-pass="dead_code_removal""

cd src
${CC} -I$CONDA_PREFIX/include -c mpi_wrapper.c
${FC} -c $FFLAGS mpi_c_bindings.f90
${FC} -c $FFLAGS mpi.f90
${FC} -c $FFLAGS psi_io.f90 --no-style-warnings --no-warnings
${FC} -c $FFLAGS --cpp --implicit-interface pot3d.F90 --no-style-warnings --no-warnings
${FC} $FFLAGS mpi_wrapper.o mpi_c_bindings.o mpi.o psi_io.o pot3d.o -o pot3d -L$CONDA_PREFIX/lib -lmpi -Wl,-rpath,$CONDA_PREFIX/lib
cp pot3d ../bin/
