#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi

printf "Testing layer decomposition: "
$1 lot_of_layers.fmm >test_output.txt

# An awk script is used to calculate the sum of the absolute values
# of the differences between the results present in the different files.
# Here the meaningful result is on the fourth column, being the result 
# of a propagation.

awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {comp=0; norm=0;} \
	FNR==NR{a[FNR]=\$4; norm+=\$4*\$4; next} \
	a[FNR]!=\$4{comp+=abs(\$4-a[FNR])} \
	END{\
		norm=sqrt(norm/FNR);\
		if (comp>1e-4*norm){ \
			exit 1 \
		}\
	}"
	
echo $awk_script > awk_scr.awk 

if awk -f awk_scr.awk test_results.txt test_reference.txt 
then
  printf " \x1b[32mOK\033[0m\n"
  diff test_results.txt test_reference.txt >test_output.txt
  rm test_results.txt test_output.txt
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail