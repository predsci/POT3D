set -ex

if [[ "$(uname)" == "Linux" ]]; then
  CC=gcc
else
  CC=clang
fi

FC=lfortran

cd src
${CC} -I$CONDA_PREFIX/include -c mpi_wrapper.c
${FC} -c mpi_c_bindings.f90
${FC} -c mpi.f90
${FC} -c psi_io.f90 --no-style-warnings --no-warnings
${FC} -c --cpp --implicit-interface pot3d.F90 --no-style-warnings --no-warnings
${FC} mpi_wrapper.o mpi_c_bindings.o mpi.o psi_io.o pot3d.o -o pot3d -L$CONDA_PREFIX/lib -lmpi -Wl,-rpath,$CONDA_PREFIX/lib
cp pot3d ../bin/
