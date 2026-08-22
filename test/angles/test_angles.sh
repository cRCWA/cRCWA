

#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to CRCWA executable.\n"
	printf "For example: $0 ../../crcwa\n"
	exit 1
fi
printf "Testing angles:              "

if [ -e "ERROR.txt" ]; then
	rm ERROR.txt
fi

$1 angles.fmm >log.txt
$1 excitation_index.fmm >>log.txt

# gnuplot checking
# plot 'test_output_index.txt' u 1:4 w p , 'test_output_index.txt' u 1:5 w l t 'ref', 'test_output_angles.txt' u 1:4 w p


if [ -e "ERROR.txt" ]; then
    printf "\033[1m\x1b[31m FAIL \033[0m\n"
    test_fail=1
else
  printf " \x1b[32mOK\033[0m\n"
  rm test_output*.txt
  rm log.txt
fi

rm Ey.txt

exit $test_fail
