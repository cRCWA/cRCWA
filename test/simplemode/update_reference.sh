#!/bin/bash

../../afmm simplemode.fmm >test_output.txt

cp test_results.txt test_reference.txt
cp test_Ex2_r_0.mode test_reference_mode_0.mode
cp test_Ey2_r_1.mode test_reference_mode_1.mode
cp test_Ex2_r_2.mode test_reference_mode_2.mode

rm test_results.txt test_output.txt test_Ex2_r_0.mode test_Ex2_r_1.mode test_Ex2_r_2.mode test_Ey2_r_0.mode test_Ey2_r_1.mode test_Ey2_r_2.mode
