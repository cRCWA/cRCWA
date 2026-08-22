#!/bin/bash

../../bin/crcwa material.fmm > log.txt

# before runing this script we should check that the interpolated
# refractive index lies on the reference curve by plotting:

# plot 'Aspnes_Si.txt' u 1:2 w l t 'ref', 'Aspnes_Si_result.txt' u 1:2 w p, 'Aspnes_Si_result.txt' u 1:3 w p
# plot 'Aspnes_Si2.txt' u 1:2 w l t 'ref', 'Aspnes_Si2_result.txt' u 1:2 w p, 'Aspnes_Si2_result.txt' u 1:3 w p
# plot 'lambda_n_k_Si.dat' u 1:2 w l t 'ref nr', 'lambda_n_k_Si.dat' u 1:3 w l t 'ref ni', 'lambda_n_k_Si_result.txt' u 1:2 w p, 'lambda_n_k_Si_result.txt' u 1:3 w p

cp Aspnes_Si_result.txt Aspnes_Si_result_reference.txt
cp Aspnes_Si2_result.txt Aspnes_Si2_result_reference.txt
cp lambda_n_k_Si_result.txt lambda_n_k_Si_result_reference.txt
cp inpstruct_result.txt inpstruct_result_reference.txt

rm log.txt

