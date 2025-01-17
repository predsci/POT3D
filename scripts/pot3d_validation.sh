#!/bin/bash

new_run=$1
valid_run=$2

cX="\033[0m"
cG="\033[32m"
cR="\033[31m"
echo="echo -e"

br1=$(grep "Energy in Br" $1)
br2=$(grep "Energy in Br" $2)

bt1=$(grep "Energy in Bt" $1)
bt2=$(grep "Energy in Bt" $2)

bp1=$(grep "Energy in Bp" $1)
bp2=$(grep "Energy in Bp" $2)

if [ "$br1" == "" ]; then
   ${echo} "Run may have ${cR}FAILED${cX} to run!"
   exit 1
fi

passed=0

if [ "$br1" != "$br2" ]; then
  passed=$(($passed + 1))
fi	

if [ "$bt1" != "$bt2" ]; then
  passed=$(($passed + 1))
fi

if [ "$bp1" != "$bp2" ]; then
  passed=$(($passed + 1))
fi

if [ $passed == 0 ]; then
  ${echo} "Run has ${cG}PASSED${cX} validation!"
else
  ${echo} "Run may have ${cR}FAILED${cX} validation!"
  echo " "
  echo "RUN1 BR: $br1"
  echo "RUN2 BR: $br2"
  echo " "
  echo "RUN1 BT: $bt1"
  echo "RUN2 BT: $bt2"
  echo " "
  echo "RUN1 BP: $bp1"
  echo "RUN2 BP: $bp2"
  # to make sure that the CI fails
  exit 1
fi

exit 0
