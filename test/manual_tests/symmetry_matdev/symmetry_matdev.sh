# testing the eigenvalues obtained for a rectangular waveguide for several 
# identical configurations.

AFMM=../../bin/afmm

# with default matrix developpement and no-symmetry setup:
echo "DEFAULT SOLVER"
echo ""
echo "default setup"
FILE=rect_wguide
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 


# with lalanne developpement = 0 and no-symmetry setup:
echo "Lalanne developpement 0"
FILE=rect_wguide_la0
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 

# with lalanne developpement = 1 and no-symmetry setup:
echo "Lalanne developpement 1"
FILE=rect_wguide_la1
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 

# with normal field developpement and no-symmetry setup:
echo "Normal field developpement"
FILE=rect_wguide_nf
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 






echo "SYMMETRY SOLVER"
echo ""
echo " without symmetry"
# with lalanne developpement = 0 and symmetry setup:
echo "Lalanne developpement 0"
FILE=rect_wguide_la0_symnn
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 

# with lalanne developpement = 1 and symmetry setup:
echo "Lalanne developpement 1"
FILE=rect_wguide_la1_symnn
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 

# with normal field developpement and symmetry setup:
echo "Normal field developpement"
FILE=rect_wguide_nf_symnn
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 

echo ""
echo " with symmetry a s"
# with lalanne developpement = 0 and symmetry setup:
echo "Lalanne developpement 0"
FILE=rect_wguide_la0_symas
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 
grep "Interesting mode" $FILE.log  > $FILE.out 

# with lalanne developpement = 1 and symmetry setup:
echo "Lalanne developpement 1"
FILE=rect_wguide_la1_symas
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 
grep "Interesting mode" $FILE.log  > $FILE.out 

# with normal field developpement and symmetry setup:
echo "Normal field developpement"
FILE=rect_wguide_nf_symas
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7
grep "Interesting mode" $FILE.log  > $FILE.out 
echo ""
echo " with symmetry s a"
# with lalanne developpement = 0 and symmetry setup:
echo "Lalanne developpement 0"
FILE=rect_wguide_la0_symsa
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7
grep "Interesting mode" $FILE.log  > $FILE.out 

# with lalanne developpement = 1 and symmetry setup:
echo "Lalanne developpement 1"
FILE=rect_wguide_la1_symsa
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 
grep "Interesting mode" $FILE.log  > $FILE.out 

# with normal field developpement and symmetry setup:
echo "Normal field developpement"
FILE=rect_wguide_nf_symsa
$AFMM $FILE.fmm > $FILE.log 
grep "Interesting mode" $FILE.log | sort -k7 
grep "Interesting mode" $FILE.log  > $FILE.out 

# REPORT:
# We can combine the mode analysis with sa and as symmetry to get the correct
# modes (depending on their symmetry).
echo ""
echo "" 
echo "combining symmetry a s and symmetry s a"
echo "some mode may not have the correct symmetry and thus will not propagate"
echo "light but they are found !!!"
# remove the modes that do not have the good symmetry
sed '1,1d' rect_wguide_la0_symsa.out >  rect_wguide_la0_symsa.out2
sed '3,3d' rect_wguide_la0_symas.out >  rect_wguide_la0_symas.out2

sed '1,1d' rect_wguide_la1_symsa.out >  rect_wguide_la1_symsa.out2
sed '3,3d' rect_wguide_la1_symas.out >  rect_wguide_la1_symas.out2

sed '2,2d' rect_wguide_nf_symsa.out >  rect_wguide_nf_symsa.out2
sed '4,4d' rect_wguide_nf_symas.out >  rect_wguide_nf_symas.out2


echo "Lalanne developpement 0"
sort -k7 rect_wguide_la0_symsa.out2 rect_wguide_la0_symas.out2
echo "Lalanne developpement 1"
sort -k7 rect_wguide_la1_symsa.out2 rect_wguide_la1_symas.out2
echo "Normal field developpement"
sort -k7 rect_wguide_nf_symsa.out2 rect_wguide_nf_symas.out2


# clean some files 
rm *.out*




