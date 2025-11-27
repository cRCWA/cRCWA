# testing all monitors for calculating the absorptance of a Si NW at 400 nm
# for several identical configurations.
# it uses the normal field developpement


AFMM=../../bin/afmm
POLYGONE=../../normal_field/polygone_c-s
REPORT=report.txt

# Creation of the structure that is a nanowire of Si studied at 500 nm.
echo "CREATING THE STRUCTURE"
echo ""
$POLYGONE nanowire.f3d 1 1 0.25 50 0.5 0.5 0 5.57 -0.387 1 0

# prepare the report:
echo "Report obtained for the calculation of absobed power in a nanowire of" > $REPORT
echo "Silicon on a non absorbing substrate" >> $REPORT
echo "" >> $REPORT
echo "filename generation generation_averaged power powerZ monitor" >> $REPORT


# with default matrix developpement and no-symmetry setup:
echo "using the normal field developpement"
echo "DEFAULT SOLVER"
echo ""
echo "default setup"
FILE=cylinder_nanowire_std
$AFMM $FILE.fmm > $FILE.log 
grep "power" $FILE.log > $FILE.out
./absorbed_power.sh $FILE $REPORT
# display the generation rate:
echo "	set pm3d map
		set terminal png size 800,800
		set cbrange [0:3e23]
		set output '$FILE.generation.png'
		splot '$FILE.generation'" > plot.plt 
gnuplot plot.plt


# with default matrix developpement and symmetry solver n n:
echo ""
echo "with symmetry n n"
FILE=cylinder_nanowire_symnn
$AFMM $FILE.fmm > $FILE.log 
grep "power" $FILE.log > $FILE.out
./absorbed_power.sh $FILE $REPORT
# display the generation rate:
echo "	set pm3d map
		set terminal png size 800,800
		set cbrange [0:3e23]
		set output '$FILE.generation.png'
		splot '$FILE.generation'" > plot.plt 
gnuplot plot.plt


# with default matrix developpement and symmetry solver a s:
echo "with symmetry a s"
FILE=cylinder_nanowire_symas
$AFMM $FILE.fmm > $FILE.log 
grep "power" $FILE.log > $FILE.out
./absorbed_power.sh $FILE $REPORT
# display the generation rate:
echo "	set pm3d map
		set terminal png size 800,800
		set cbrange [0:3e23]
		set output '$FILE.generation.png'
		splot '$FILE.generation'" > plot.plt 
gnuplot plot.plt





# with default matrix developpement and symmetry solver s a:
#echo "with symmetry s a"


# printing the report
echo ""
echo "REPORT"
cat $REPORT

# clean useless files
rm *.out
rm plot.plt

