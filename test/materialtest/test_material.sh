#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing material:            "

$1 material.fmm >log.txt

# gnuplot checking
# plot 'Aspnes_Si.txt' u 1:2 w l t 'ref', 'Aspnes_Si_result.txt' u 1:2 w p, 'Aspnes_Si_result.txt' u 1:3 w p
# plot 'Aspnes_Si2.txt' u 1:2 w l t 'ref', 'Aspnes_Si2_result.txt' u 1:2 w p, 'Aspnes_Si2_result.txt' u 1:3 w p
# plot 'lambda_n_k_Si.dat' u 1:2 w l t 'ref nr', 'lambda_n_k_Si.dat' u 1:3 w l t 'ref ni', 'lambda_n_k_Si_result.txt' u 1:2 w p, 'lambda_n_k_Si_result.txt' u 1:3 w p


if diff Aspnes_Si_result.txt Aspnes_Si_result_reference.txt >test_output.txt && diff Aspnes_Si2_result.txt Aspnes_Si2_result_reference.txt >>test_output.txt && diff lambda_n_k_Si_result.txt lambda_n_k_Si_result_reference.txt && diff inpstruct_result.txt inpstruct_result_reference.txt  >>test_output.txt
then
  printf " \x1b[32mOK\033[0m\n"
  rm log.txt test_output.txt
else
  test_fail=1
  printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail


