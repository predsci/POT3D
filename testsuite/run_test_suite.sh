#!/bin/bash

trap ctrl_c INT

function ctrl_c() {
  ${echo} "${cR}==> Caught CTRL-C, shutting down!${cX}"
  exit 1
}

function display_help {
echo " _____   ____ _______ ____  _____    "
echo "|  __ \ / __ \__   __|___ \|  __ \   "
echo "| |__) | |  | | | |    __) | |  | |  "
echo "|  ___/| |  | | | |   |__ <| |  | |  "
echo "| |    | |__| | | |   ___) | |__| |  "
echo "|_|     \____/  |_|  |____/|_____/   "
echo " "
echo "TEST SUITE v1.0.0"
echo "USAGE:   ./run_test_suite.sh"
echo ""
echo "By default, the above command will run the test suite on all"
echo "default tests using the 'pot3d' executable from the '../bin/'"
echo "folder (assuming it has been built)."
echo ""
echo "OPTIONS:"
echo ""
echo "Non-flag options '-opt=' are specified as '-opt=<OPTION>'."
echo ""
echo "-nochecksetup   Don't check the environment."
echo ""
echo "-pot3dexe=      Use this to run the testsuite on a specific pot3d executable."
echo "                This should be a full path and is useful for development tests."
echo ""
echo "-mpicall=       Use this to run the testsuite with a specific mpi calling mechanism."
echo "                The call should end with '-np' or equivalent as the # of ranks will"
echo "                be placed after the call. "
echo "                For example, for GCC+OpenMPI, one could use: mpicall='mpirun -bind-to socket -np'"
echo "                The default is 'mpirun -np'."
echo ""
echo "-np=            Number of MPI ranks to use."
echo ""
echo "-test=          Comma-seperated list of subset of tests to run."
echo "                Also can be used to run non-standard/experimental tests."
echo ""
echo "-nocleanup      By default, only the initial and final conditions of a run are "
echo "                kept in the run folders. Set this to keep the full run."
echo ""
echo "-cleanup        Remove all previous test suite runs (if they were run with -nocleanup)"
echo ""
echo "-norun          Does not run the POT3D code.  Checks for a previous run and compares if found."
echo ""
echo "-nocompare      Do not compare the runs to their reference runs."
echo ""
echo "-compareprec=   Set the precision for the comparisons (decimal place)"
echo "                Default is 5."
echo ""
echo "-nocolor        Disable color text output."
echo ""
}

#Set number of processors to use for testsuite:
#(All tests use the same number of procs for now).
np=1
norun=0
nocompare=0
compareprec=5
novis=0
nocleanup=0
cleanup=0
nochecksetup=0
nocolor=0
setrefdata=0
pc=1
pot3dexe="pot3d"
mpicall="mpirun -np"

AVAIL_TEST_RUNS_LIST="
potential_field_source_surface
potential_field_current_sheet
open_field
"

TEST_RUNS_LIST=${AVAIL_TEST_RUNS_LIST}

for i in "$@"
do
case $i in
    -norun)
    norun=1
    ;;
    -nocompare)
    nocompare=1
    ;;
    -compareprec=*)
    compareprec="${i#*=}"
    ;;    
    -nocleanup)
    nocleanup=1
    ;;
    -nochecksetup)
    nochecksetup=1
    ;;
    -nocolor)
    nocolor=1
    ;;
    -setrefdata)
    setrefdata=1
    ;;
    -test=*)
    TEST_RUNS_LIST="${i#*=}"
    TEST_RUNS_LIST="${TEST_RUNS_LIST//','/' '}"
    ;;
    -pot3dexe=*)
    pot3dexe="${i#*=}"
    ;;
    -mpicall=*)
    mpicall="${i#*=}"
    ;;
    -np=*)
    np="${i#*=}"
    ;;
    -pc=*)
    pc="${i#*=}"
    ;;
    -cleanup)
    cleanup=1
    norun=1
    nocompare=1
    nochecksetup=1
    ;;
    -h)
    display_help
    exit 0
    ;;
    --help)
    display_help
    exit 0
    ;;    
    *)
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "ERROR!  Unknown option: $i"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    display_help
    exit 1    
    ;;
esac
done

# ****** Get test suite parameters ******
WD=${PWD}
BINDIR=${WD}/../bin
SRCDIR=${WD}/../src
ROOTDIR=${WD}/..
TSLOG=${RESULTSDIR}/testsuite.log

if [ ${nocolor} == 0 ]
then
  cX="\033[0m"
  cR="\033[1;31m"
  cB="\033[1;34m"
  cG="\033[32m"
  cC="\033[1;96m"
  cM="\033[35m"
  cY="\033[1;93m"
  Bl="\033[1;5;96m"
  echo="echo -e"
else
  cX=
  cR=
  cB=
  cG=
  cC=
  cM=
  cY=
  Bl=
  echo="echo"
fi

${echo} " "
${echo} "${cB}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
${echo} "${cY}   _____   ____ _______ ____  _____                                   ${cX}"
${echo} "${cY}  |  __ \ / __ \__   __|___ \|  __ \                                  ${cX}"
${echo} "${cY}  | |__) | |  | | | |    __) | |  | |                                 ${cX}"
${echo} "${cY}  |  ___/| |  | | | |   |__ <| |  | |                                 ${cX}"
${echo} "${cY}  | |    | |__| | | |   ___) | |__| |                                 ${cX}"
${echo} "${cY}  |_|     \____/  |_|  |____/|_____/                                  ${cX}"
${echo} "${cY}  ___ ____ ____ ___    ____ _  _ _ ___ ____                           ${cX}"
${echo} "${cY}   |  |___ [__   |     [__  |  | |  |  |___                           ${cX}"
${echo} "${cY}   |  |___ ___]  |     ___] |__| |  |  |___                           ${cX}"
${echo} "${vY}                                                                      ${cX}"
${echo} "${cB}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
${echo} "Welcome to the POT3D test suite!"
#
# ****** Test for correct prerequisites and environment ******
#
${echo} "Checking file structure..."

#
# ****** Test for correct prerequisites and environment ******
#
if [ ${nochecksetup} == 0 ]
then
#
# Check that the POT3D bin directory is in the user's path, if not, add it.
#
  ${echo} "==> Checking PATH...."
  PTEST=$(which pot3d_validation.sh)
  if [ -z "${PTEST}" ]
  then
    ${echo} "${cY}==> WARNING: POT3D bin is not in the PATH!${cX}"
    ${echo} "${cY}==> Appending ${BINDIR} to PATH...${cX}"
    PATH="${BINDIR}:${PATH}"
  fi
  PTEST=$(which pot3d)
  if [ -z "${PTEST}" ]; then
    ${echo} "${cR}==> ERROR! POT3D not in bin, perhaps it has not been built?${cX}"
    exit 1
  fi
  ${echo} "${cG}==> Everything seems OK to run POT3D test suite!${cX}"
fi
#
# ****** Check that user is running the tests and really wants to set ref data:
#
if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 1 ]
then
  ${echo} "${cR} ==> ERROR!  You are trying to set reference data without running the tests!${cX}"
  exit 1
fi
if [ ${setrefdata} -eq 1 ]
then
  ${echo} "${cR}╔═╗┌─┐┌┬┐  ╔╗╔┌─┐┬ ┬  ╦═╗┌─┐┌─┐┌─┐┬─┐┌─┐┌┐┌┌─┐┌─┐  ╔╦╗┌─┐┌┬┐┌─┐${cX}"
  ${echo} "${cR}╚═╗├┤  │   ║║║├┤ │││  ╠╦╝├┤ ├┤ ├┤ ├┬┘├┤ ││││  ├┤    ║║├─┤ │ ├─┤${cX}"
  ${echo} "${cR}╚═╝└─┘ ┴   ╝╚╝└─┘└┴┘  ╩╚═└─┘└  └─┘┴└─└─┘┘└┘└─┘└─┘  ═╩╝┴ ┴ ┴ ┴ ┴${cX}"
  read -p "==> Setting new reference data after runs complete...are you SURE? (y/n):" yn
  if [ ${yn} = "y" ]
  then
    read -p "==> Are you really really really SURE?! (y/n):" yn
    if [ ${yn} = "n" ]
    then
      ${echo} "==> ${cR}Aborting!${cX}"
      exit 0
    fi
  else
    ${echo} "==> ${cR}Aborting!${cX}"
    exit 0
  fi
fi

########################################################################
########################################################################
##
## ****** Start loop through test problems ******
##
########################################################################
########################################################################

Ti=0
for TESTNAME in ${TEST_RUNS_LIST}
do
  Ti=$((${Ti}+1))
#
# ****** Make sure test is in the test suite:
#
  ${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
  testok=0
  for run_test in ${AVAIL_TEST_RUNS_LIST}; do
    if [[ ${run_test} == ${TESTNAME} ]]; then
      testok=1
      break
    fi
  done
  if [ ${testok} -eq 0 ]
  then
    ${echo} "${cR}==> TEST ${cX}${cM}${TESTNAME}${cX}${cR} is not a valid test in the testsuite!${cX}"
    continue
  fi
  ${echo} "STARTING TEST ${cM}${TESTNAME}${cX}"
  ${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
  ${echo} "==> Gathering test information..."
# ****** Get directories:
  RUNDIR=${WD}/${TESTNAME}/run
  REFDIR=${WD}/${TESTNAME}/reference
  INPUTDIR=${WD}/${TESTNAME}/input

  if [ ! -d ${RUNDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Run directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -d ${REFDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Reference directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -d ${INPUTDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Input directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -e ${INPUTDIR}/pot3d.dat ]
  then
    ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} does not have an input file!${cX}"
    exit 1
  fi

  cd ${RUNDIR}
  
  if [ ${cleanup} == 1 ]
  then
    ${echo} "==> Clearing run directory..."
    rm -fr ${RUNDIR}/*
  fi

  if [ ${norun} == 0 ]
  then

    if [ -e ${RUNDIR}/pot3d.out ]
    then
      ${echo} "==> Clearing run directory..."
      rm -fr ${RUNDIR}/*
    fi

    ${echo} "==> Copying input files..."
    cp ${INPUTDIR}/* ${RUNDIR}/ 2>/dev/null

    cd ${RUNDIR}

    # Set preconditioner:
    (sed -i '' "s/\([[:space:]]*\)ifprec=1\([[:space:]]*\)/\1ifprec=$pc\2/g" pot3d.dat 2>/dev/null \
    || sed -i "s/\([[:space:]]*\)ifprec=1\([[:space:]]*\)/\1ifprec=$pc\2/g" pot3d.dat)

    ${echo} "======================================================="
    ${echo} "${cB}==> RUNNING POT3D${cX}"
    ${echo} "======================================================="
    ${echo} "==> Running pot3d with command:"
    ${echo} "==> ${mpicall} $np ${pot3dexe} pot3d.dat 1>pot3d.log 2>pot3d.err"
    ${mpicall} $np ${pot3dexe} pot3d.dat 1>pot3d.log 2>pot3d.err

    # Check that a completed run exists in the run folder
    if [ ! -e ${RUNDIR}/br.h5 ]
    then
      if [ ${norun} == 0 ]
      then
        ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} did not seem to run correctly!${cX}"
        ${echo} "${cR}!!!> pot3d.log contents: ${cX}"
        cat ${RUNDIR}/pot3d.log
        ${echo} "${cR}!!!> pot3d.err contents: ${cX}"
        cat ${RUNDIR}/pot3d.err
        ${echo} "${cR}!!!> Check the run folder: ${RUNDIR} ${cX}"
        exit 1
      fi
    fi

  # Get timing data:
    TIME_RUN_TMP=($(grep "Average time used per proc" ${RUNDIR}/timing.out))
    TIME_RUN_TMP=${TIME_RUN_TMP[6]}
    TIME_RUN[${Ti}]=$(printf "%8.3f" ${TIME_RUN_TMP})
    ${echo} "${cG}==> Test completed in:               ${TIME_RUN[${Ti}]} seconds.${cX}"
    
    if [ ${nocompare} == 0 ]
    then
      TIME_REF_TMP=($(grep "Average time used per proc" ${REFDIR}/timing.out))
      TIME_REF_TMP=${TIME_REF_TMP[6]}

      SPEEDUP_TMP=`python3 -c "print(${TIME_REF_TMP}/${TIME_RUN_TMP})"`
      TIME_REF[${Ti}]=$(printf "%8.3f" ${TIME_REF_TMP})
      SPEEDUP[${Ti}]=$(printf "%5.2f" ${SPEEDUP_TMP})

      ${echo} "${cG}==> Speedup compared to reference run: ${SPEEDUP[${Ti}]} x${cX}"
    fi
  fi

  if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 0 ]
  then
    ${echo} "${cR}=======================================================${cX}"
    ${echo} "${cR}==> SETTING REFERENCE DATA FOR RUN${cX}"
    ${echo} "${cR}=======================================================${cX}"

    ${echo} "${cR}==> Removing old reference data...${cX}"
    rm -fr ${REFDIR}/* 2>/dev/null
    ${echo} "${cR}==> Copying current run data into reference directory...${cX}"
    cp ${RUNDIR}/* ${REFDIR}/
    rm -f ${REFDIR}/*.h5
  fi

#
# ****** Compare run data.
#
  if [ ${nocompare} == 0 ]
  then
    ${echo} "======================================================="
    ${echo} "${cB}==> COMPARING RUN DATA TO REFERENCE DATA${cX}"
    ${echo} "======================================================="
    ${echo} "==> Running comparison..."
    ${echo} -n "${cR}"
    pot3d_validation.sh ${RUNDIR}/pot3d.out ${REFDIR}/pot3d.out
    PASS_FAIL[${Ti}]=$?
    ${echo} -n "${cX}"

    if [ "${PASS_FAIL[${Ti}]}" = "0" ]
    then
      ${echo} "${cG}==> Test seems to have PASSED!${cX}"
    else
      ${echo} "${cR}==> Test seems to have FAILED!${cX}"
      ${echo} "${cR}==> ${PASS_FAIL[${Ti}]} ${cX}"
      PASS_FAIL[${Ti}]=1
    fi
  fi
#
# ****** Cleanup run data.
#
  if [ ${nocleanup} == 0 ]
  then
    if [ ${norun} == 0 ]
    then
      ${echo} "======================================================="
      ${echo} "${cB}==> CLEANING RUN DATA${cX}"
      ${echo} "======================================================="
      ${echo} "==> Removing files from run data..."
      rm -fr ${RUNDIR}/*
    fi
  fi
done
${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"


# Display summary and timing results.
if [ ${nocompare} == 0 ]
then
  ${echo} "${cY}===========================================================================${cX}"
  ${echo} "${cY}Summary of test results:${cX}"
  ${echo} "${cY}===========================================================================${cX}"

  Ti=0
    ${echo} "$(printf "%-35s  %9s  %8s  %8s  %7s" "Test name" "PASS/FAIL" "Run-time" "Ref-time" "Speedup")"
  ${echo} "${cY}===========================================================================${cX}"    
  for TESTNAME in ${TEST_RUNS_LIST}
  do
    Ti=$((${Ti}+1))
    passfailcomp=( ${PASS_FAIL[${Ti}]} )
    if [ "${passfailcomp[0]}" = "1" ]; then
      pf="${cR}FAIL     ${cX}"
    else
      pf="${cG}PASS     ${cX}"
    fi
    ${echo} "$(printf "%-35s  %9s  %8s  %8s  %7s" "${TESTNAME}" "${pf}" "${TIME_RUN[${Ti}]}" "${TIME_REF[${Ti}]}" "${SPEEDUP[${Ti}]}")"
  done
fi

${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"

exit 0



