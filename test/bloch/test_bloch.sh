#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing Bloch mode calc. TE: "

$1 bloch_TE.fmm >test_output.txt
mv propagEx.txt test_results_TE.txt


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
		if (comp>2e-4*norm){ \
			exit 1 \
		}\
	}"
	
echo $awk_script > awk_scr.awk 

if awk -f awk_scr.awk test_results_TE.txt test_reference_TE.txt
then
  printf " \x1b[32mOK\033[0m\n"
  rm test_results_TE.txt test_output.txt
else
  test_fail=1
  diff test_results_TE.txt test_reference_TE.txt >test_output.txt
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

printf "Testing Bloch mode calc. TM: "

$1 bloch_TM.fmm >test_output.txt
mv propagEy.txt test_results_TM.txt

if awk -f awk_scr.awk test_results_TM.txt test_reference_TM.txt
then 
  printf " \x1b[32mOK\033[0m\n"
  rm test_results_TM.txt test_output.txt
else
  test_fail=1
  diff test_results_TM.txt test_reference_TM.txt >test_output.txt
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail
