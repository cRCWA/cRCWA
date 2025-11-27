# testing the representation of guided modes obtained for a cylindrical waveguide
# for several identical configurations.
# it uses the normal field developpement

AFMM=../../bin/afmm
POLYGONE=../../normal_field/polygone_c-s
# number of mode that we will study
NBMODES=6	# multiple of 2

# Creation of the structure that is a nanowire of Si studied at 500 nm.
echo "CREATING THE STRUCTURE"
echo ""
$POLYGONE wguide.f3d 1 1 0.25 50 0.5 0.5 0 4.297 -0.0729 1 0

# with default matrix developpement and no-symmetry setup:
echo "using the normal field developpement"
echo "DEFAULT SOLVER"
echo ""
echo "default setup"
FILE=cylinder_wguide_std
$AFMM $FILE.fmm > $FILE.log 
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=,
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=, > $FILE.out
grep "Interesting mode" $FILE.log | sort -k7 -r
grep "Interesting mode" $FILE.log | sort -k7 -r > $FILE.out  

#save pictures:
./savemode.sh $NBMODES $FILE
# clean the mode files
rm *.mode

# with default matrix developpement and symmetry solver n n:
echo ""
echo "with symmetry n n"
FILE=cylinder_wguide_symnn
$AFMM $FILE.fmm > $FILE.log 
# sort from the smaller imaginary refrative intex to the larger
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=,
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=, > $FILE.out 
grep "Interesting mode" $FILE.log | sort -k7 -r
grep "Interesting mode" $FILE.log | sort -k7 -r > $FILE.out  
#save pictures
./savemode.sh $NBMODES $FILE
# clean the mode files
rm *.mode


# with default matrix developpement and symmetry solver a s:
echo ""
echo "with symmetry a s"
FILE=cylinder_wguide_symas
$AFMM $FILE.fmm > $FILE.log 
# sort from the smaller imaginary refrative intex to the larger
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=,
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=, > $FILE.out 
grep "Interesting mode" $FILE.log | sort -k7 -r
grep "Interesting mode" $FILE.log | sort -k7 -r > $FILE.out  
#save pictures
nbmodes=`expr $NBMODES / 2`
echo $nbmodes
./savemode.sh `expr $NBMODES / 2` $FILE
# clean the mode files
rm *.mode



# with default matrix developpement and symmetry solver s a:
echo ""
echo "with symmetry s a"
FILE=cylinder_wguide_symsa
$AFMM $FILE.fmm > $FILE.log 
# sort from the smaller imaginary refrative intex to the larger
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=,
#grep "Interesting mode" $FILE.log | sort -k3 --field-separator=, > $FILE.out 
grep "Interesting mode" $FILE.log | sort -k7 -r
grep "Interesting mode" $FILE.log | sort -k7 -r > $FILE.out  
#save pictures
./savemode.sh `expr $NBMODES / 2` $FILE
# clean the mode files
rm *.mode



# clean useless files
rm *.out

