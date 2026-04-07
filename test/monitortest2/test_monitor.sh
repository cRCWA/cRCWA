#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing monitor 2:           "

$1 monitor.fmm >log.txt

# printf "result output: "
# cat test_output.txt

# An awk script is used to calculate the sum of the absolute values
# of the differences between the results present in the different files.
# Here the meaningful result is on the third column, being the result 
# of a mode calculation (real part).
# awk 'NR==FNR{a[$1]=$2;next} sqrt((a[$1]-$3)^2) && sqrt((a[$1]-$4)^2) {print $1":"$2":"a[$1]}' test_output.txt test_output_reference.txt


awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {error=0;} \
	FNR==NR{id[FNR]=\$1; p[FNR]=\$2; pz[FNR]=\$3; m1[FNR]=\$4; m2[FNR]=\$5; m3[FNR]=\$6; next} \
	abs(p[\$1]-\$2)>1e-18 || abs(pz[\$1]-\$3)>1e-18 || abs(m1[\$1]-\$4)>1e-18 || abs(m2[\$1]-\$5)>1e-18 || abs(m3[\$1]-\$6)>1e-18 \
	{error=1; print \$1\" \" p[\$1] \" \" pz[\$1] \" \" m1[\$1]\" \"m2[\$1] \" \" m3[\$1] \" vs \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$6} \
	END{exit error }"
	

echo $awk_script > awk_scr.awk 

if awk -f awk_scr.awk test_output.txt test_output_reference.txt 
then
  printf " \x1b[32mOK\033[0m\n"
  rm Ex*.txt
  rm Ey*.txt
  rm Hx*.txt
  rm Hy*.txt
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail
