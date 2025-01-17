#!/bin/bash

POT3D_HOME=$PWD
TEST="validation"

cp ${POT3D_HOME}/testsuite/${TEST}/input/* ${POT3D_HOME}/testsuite/${TEST}/run/
cd ${POT3D_HOME}/testsuite/${TEST}/run

echo "Running POT3D with 1 MPI rank..."
mpiexec -np 1 ${POT3D_HOME}/bin/pot3d 1> pot3d.log 2>pot3d.err
echo "Done!"

runtime=($(tail -n 5 timing.out | head -n 1))
echo "Wall clock time:                ${runtime[6]} seconds"
echo " "

#Validate run:
${POT3D_HOME}/scripts/pot3d_validation.sh pot3d.out ${POT3D_HOME}/testsuite/${TEST}/validation/pot3d.out
if [ $? -ne 0 ]; then
  echo "Validation failed for 1 MPI rank. Exiting..."
  exit 1
fi

rm pot3d.log pot3d.out timing.out
echo " "
echo "Running POT3D with 2 MPI ranks..."
mpiexec -np 2 ${POT3D_HOME}/bin/pot3d 1> pot3d.log 2>pot3d.err
echo "Done!"

runtime=($(tail -n 5 timing.out | head -n 1))
echo "Wall clock time:                ${runtime[6]} seconds"
echo " "

#Validate run:
${POT3D_HOME}/scripts/pot3d_validation.sh pot3d.out ${POT3D_HOME}/testsuite/${TEST}/validation/pot3d.out
if [ $? -ne 0 ]; then
  echo "Validation failed for 1 MPI rank. Exiting..."
  exit 1
fi
echo " "

