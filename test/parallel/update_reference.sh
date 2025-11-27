#!/bin/bash

../../afmm lot_of_layers.fmm >test_output.txt

cp test_results.txt test_reference.txt

rm test_results.txt