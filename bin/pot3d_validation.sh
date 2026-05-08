#!/bin/bash
new_run=$1
valid_run=$2
decimals=${3:-7}

extract_num() {
    echo "$1" | grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' | tail -1
}

round() {
    awk -v n="$decimals" -v val="$1" '
        BEGIN {
            factor = 10^n
            printf "%.0f\n", val * factor + 0.5 * (val > 0 ? 1 : -1)
        }' | awk -v n="$decimals" '{ printf "%.*f\n", n, $1 / (10^n) }'
}

br1=$(grep "Energy in Br" "$new_run" | head -1)
br2=$(grep "Energy in Br" "$valid_run" | head -1)
bt1=$(grep "Energy in Bt" "$new_run" | head -1)
bt2=$(grep "Energy in Bt" "$valid_run" | head -1)
bp1=$(grep "Energy in Bp" "$new_run" | head -1)
bp2=$(grep "Energy in Bp" "$valid_run" | head -1)

if [ "$br1" == "" ]; then
   exit 1
fi

num_br1=$(extract_num "$br1")
num_br2=$(extract_num "$br2")
num_bt1=$(extract_num "$bt1")
num_bt2=$(extract_num "$bt2")
num_bp1=$(extract_num "$bp1")
num_bp2=$(extract_num "$bp2")

passed=0

if [ "$(round "$num_br1")" != "$(round "$num_br2")" ]; then passed=$((passed + 1)); fi
if [ "$(round "$num_bt1")" != "$(round "$num_bt2")" ]; then passed=$((passed + 1)); fi
if [ "$(round "$num_bp1")" != "$(round "$num_bp2")" ]; then passed=$((passed + 1)); fi

if [ $passed -eq 0 ]; then
    exit 0
else
    echo "Test may have FAILED (compared to $decimals decimal places):"
    echo " "
    echo "RUN1 BR: $br1  =>  $num_br1"
    echo "RUN2 BR: $br2  =>  $num_br2"
    echo " "
    echo "RUN1 BT: $bt1  =>  $num_bt1"
    echo "RUN2 BT: $bt2  =>  $num_bt2"
    echo " "
    echo "RUN1 BP: $bp1  =>  $num_bp1"
    echo "RUN2 BP: $bp2  =>  $num_bp2"
    exit 1
fi
