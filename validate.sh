#!/bin/bash

if [ -z "$1" ]; then
    NP=1
else
    NP="$1"
fi

POT3D_HOME=$PWD
TEST="validation"

cp ${POT3D_HOME}/testsuite/${TEST}/input/* ${POT3D_HOME}/testsuite/${TEST}/run/
cd ${POT3D_HOME}/testsuite/${TEST}/run

echo "Running POT3D with $NP MPI rank..."
mpiexec -np $NP ${POT3D_HOME}/bin/pot3d 1> pot3d.log 2>pot3d.err
echo "Done!"

runtime=($(tail -n 5 timing.out | head -n 1))
echo "Wall clock time:                ${runtime[6]} seconds"
echo " "

#Validate run:
${POT3D_HOME}/bin/pot3d_validation.sh pot3d.out ${POT3D_HOME}/testsuite/${TEST}/validation/pot3d.out

echo " "

