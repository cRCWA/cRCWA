#!/bin/bash

../../afmm bloch_TE.fmm >test_output.txt
mv propagEx.txt test_reference_TE.txt

../../afmm bloch_TM.fmm >test_output.txt
mv propagEy.txt test_reference_TM.txt
