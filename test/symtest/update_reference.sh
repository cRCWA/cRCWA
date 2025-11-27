#!/bin/bash

../../afmm pyramid.fmm >test_output.txt
mv test_results_pr.txt test_reference_pr.txt
mv test_results_AR.txt test_reference_AR.txt

