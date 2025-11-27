#!/bin/bash

../../afmm testorientation.fmm >test_output.txt

mv refractive_index.txt refractive_index_ref.txt 
mv refractive_index_im.txt refractive_index_im_ref.txt 
mv propag.txt propag_ref.txt
