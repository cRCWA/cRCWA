#!/bin/bash

../../afmm calc.fmm >test_output.txt

cp test_results.txt test_reference.txt

rm test_results.txt