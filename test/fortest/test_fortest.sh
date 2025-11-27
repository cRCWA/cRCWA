#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi


printf "Testing for command (2 loops)"
$1 fortest.fmm >test_output.txt

if diff test_results.txt test_reference.txt >test_output.txt
then
  printf " \x1b[32mOK\033[0m\n"
  rm test_results.txt test_output.txt
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail
