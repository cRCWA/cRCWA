#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing orientation of fft:  "

$1 testorientation.fmm >test_output.txt

# An awk script is used to calculate the sum of the absolute values
# of the differences between the results present in the different files.
# Here the meaningful result is on the third column, being the result 
# of a mode calculation (real part).

awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {comp=0; norm=0;} \
	FNR==NR{a[FNR]=\$3; norm+=\$3*\$3; next} \
	a[FNR]!=\$3{comp+=abs(\$3-a[FNR])} \
	END{\
		norm=sqrt(norm/FNR);\
		if (comp>1e-6*norm){ \
			exit 1 \
		}\
	}"
	
echo $awk_script > awk_scr.awk 


if awk -f awk_scr.awk refractive_index.txt refractive_index_ref.txt && awk -f awk_scr.awk refractive_index_im.txt refractive_index_im_ref.txt && awk -f awk_scr.awk propag.txt propag_ref.txt
then
  printf " \x1b[32mOK\033[0m\n"
  rm refractive_index.txt refractive_index_im.txt propag.txt
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail