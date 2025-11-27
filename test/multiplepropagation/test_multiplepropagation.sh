#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi

printf "Testing multi-propagation:   "
$1 multiplepropagation.fmm >test_output.txt

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

awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {comp=0; norm=0;} \
	NR>2 && FNR==NR{a[FNR]=\$3; norm+=\$3*\$3; next} \
	a[FNR]!=\$3{comp+=abs(\$3-a[FNR])} \
	END{\
		norm=sqrt(norm/FNR);\
		if (comp>1e-4*norm){ \
			exit 1 \
		}\
	}"
	
echo $awk_script > awk_scr2.awk 

if awk -f awk_scr.awk test_results_x.txt test_reference_x.txt \
	&& awk -f awk_scr.awk test_results_z.txt test_reference_z.txt  \
	&&  awk -f awk_scr.awk test_results_y.txt test_reference_y.txt 
then
  printf " \x1b[32mOK\033[0m\n"
  rm test_results_x.txt test_results_y.txt test_results_z.txt test_output.txt
  rm modes_2_Ey2_m_0.mode modes_2_Ey2_m_1.mode modes_2_Ey2_m_2.mode
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail