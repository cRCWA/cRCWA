#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing simple mode search:  "

$1 simplemode.fmm >test_output.txt

# An awk script is used to calculate the sum of the absolute values
# of the differences between the results present in the different files.
# Here the meaningful result is on the third column, being the result 
# of a mode calculation (real part).

awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {comp=0; norm=0;} \
	FNR>1&&FNR==NR{a[FNR]=\$3; norm+=\$3*\$3; next} \
	FNR>1&&a[FNR]!=\$3{comp+=abs(\$3-a[FNR])} \
	END{\
		norm=sqrt(norm/FNR);\
		if (comp>1e-6*norm){ \
			exit 1 \
		}\
	}"
	
echo $awk_script > awk_scr.awk 


if diff test_results.txt test_reference.txt >test_output.txt && awk -f awk_scr.awk test_Ex2_r_0.mode test_reference_mode_0.mode && awk -f awk_scr.awk test_Ey2_r_1.mode test_reference_mode_1.mode
then
  printf " \x1b[32mOK\033[0m\n"
  rm test_results.txt test_output.txt test_Ex2_r_0.mode test_Ex2_r_1.mode test_Ex2_r_2.mode test_Ey2_r_0.mode test_Ey2_r_1.mode test_Ey2_r_2.mode
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail