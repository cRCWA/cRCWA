#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing pl. wave exc. angles "
                                     
$1 1D_test_angle.fmm >test_output.txt

if diff test_results.txt test_reference.txt >>test_output.txt
then
  printf " \x1b[32mOK\033[0m\n"
  rm test_results.txt test_output.txt resu
else
  test_fail=1
  diff test_results.txt test_reference.txt >>test_output.txt
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail