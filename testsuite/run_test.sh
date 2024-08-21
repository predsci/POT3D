#!/bin/bash

if [ -z "$1" ]
then
    echo "Warning!  You should supply the test name and the number of MPI ranks as arguments:"
    echo "Example:  ./run_test validation 1"
    echo ""
    echo "==> Setting test name to 'validation'"
    echo ""
    TEST="validation"
else
    TEST="$1"
fi

if [ -z "$2" ]
then
    echo "Warning!  You should supply the test name and the number of MPI ranks as arguments:"
    echo "Example:  ./run_test validation 1" 
    echo "" 
    echo "==> Setting number of MPI ranks to 1"
    echo ""
    NP="1"
else
    NP="$2"
fi

POT3D_HOME=${PWD}/..

cp ${POT3D_HOME}/testsuite/${TEST}/input/* ${POT3D_HOME}/testsuite/${TEST}/run/
cd ${POT3D_HOME}/testsuite/${TEST}/run

echo "Running POT3D for test $TEST with $NP MPI rank(s)..."
mpiexec -np $NP ${POT3D_HOME}/bin/pot3d > pot3d.log
echo "Done!"

runtime=($(tail -n 5 timing.out | head -n 1))
echo "Wall clock time:                ${runtime[6]} seconds"
echo " "

#Validate run:
${POT3D_HOME}/scripts/pot3d_validation.sh pot3d.out ${POT3D_HOME}/testsuite/${TEST}/validation/pot3d.out

echo " "

