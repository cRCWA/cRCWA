#!/bin/bash

../../afmm calculatedz.fmm >test_output.txt

cp test_results_x.txt test_reference_x.txt
cp test_results_y.txt test_reference_y.txt
cp test_results_z.txt test_reference_z.txt


rm test_results_x.txt test_results_y.txt test_results_z.txt