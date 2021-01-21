#!/bin/bash

POT3D_HOME=${PWD}/..
TEST="$1"
NP="$2"

#Check for test existing...

cp ${POT3D_HOME}/testsuite/${TEST}/input/* ${POT3D_HOME}/testsuite/${TEST}/run/
cd ${POT3D_HOME}/testsuite/${TEST}/run

echo "Running POT3D with $NP MPI rank..."
mpiexec -np $NP ${POT3D_HOME}/bin/pot3d > pot3d.log
echo "Done!"

runtime=($(tail -n 5 timing.out | head -n 1))
echo "Wall clock time:                ${runtime[6]} seconds"
echo " "

#Validate run:
${POT3D_HOME}/scripts/pot3d_validation.sh pot3d.out ${POT3D_HOME}/testsuite/${TEST}/validation/pot3d.out

echo " "

