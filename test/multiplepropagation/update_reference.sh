#!/bin/bash

../../afmm multiplepropagation.fmm >test_output.txt

cp test_results_x.txt test_reference_x.txt
cp test_results_y.txt test_reference_y.txt
cp test_results_z.txt test_reference_z.txt

mv modes_2_Ey2_m_0.mode modes_2_Ey2_m_0.mode.reference
mv modes_2_Ey2_m_1.mode modes_2_Ey2_m_1.mode.reference
mv modes_2_Ey2_m_2.mode modes_2_Ey2_m_2.mode.reference



rm test_results_x.txt test_results_y.txt test_results_z.txt