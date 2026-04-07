

#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to CRCWA executable.\n"
	printf "For example: $0 ../../crcwa\n"
	exit 1
fi
printf "Testing drawing:              "

draw="../../bin/draw"
o2g="../../bin/o2g"
start_time=$(date +%s%3N) 



# generate squares with radial normal field with centered origin
$draw --legacy square1.rid 2 1 0.3535 4 0.5 0.5 45 1 0 2 0 > log.txt
$draw --legacy -p square2.rid 2 1 4 .75 .75 .25 0.75 0.25 0.25 0.75 0.25 1 0 2 0 >> log.txt
$draw --legacy -pr square3.rid 2 1 4 .75 0.25 0.25 0.25 0.25 0.75 0.75 0.75 1 0 2 0 >> log.txt

$draw square4.rid 1 2 1 p 1 0 4 .25 -.25 -.25 -0.25 -0.25 -0.75 0.25 -0.75 2 0  >> log.txt
$draw square5.rid 1 2 1 c 1 0 0.3535 4 0 -0.5 0 2 0  >> log.txt
$draw square6.rid 1 2 1 r 1 0 .5 .5 0 -0.5 2 0  >> log.txt

# generate rectangle and triangle with radial normal field in diagonal direction
# to one of the rectangle edge
$draw --legacy -pr -o .7 .7 rectangle1.rid 1.5 2 3 0.25 0.5 0.5 0.25 0.5 0.75 2 0 4 .9 0.25 0.2 0.25 0.2 0.75 0.9 0.75 3 0 1 0 >> log.txt
$draw rectangle2.rid 1 1.5 2 p 2 0 3 -0.25 -0.25 0 -0.5 0 0 p 3 0 4 .4 -0.5 -0.3 -0.5 -0.3 0 0.4 0 1 0 -o .7 .7 >> log.txt

# generate circle and hexagon 
$draw --legacy circle_2vortexes1.rid 1 2 0.15 50 0.5 0.25 0 2 0 0.2 8 0.5 0.75 0 3 0 1 0 >> log.txt
$draw circle_2vortexes2.rid 1 1 2 c 2 0 0.15 50 0 -0.25 0 c 3 0 0.2 8 0 0.25 22.5 1 0 >> log.txt

# generate circle and hexagone with 1 vortex (origin of vector take inside hexagone)
$draw --legacy -c 0.5 0.8 circle_1vortex1.rid 1 2 0.15 50 0.5 0.25 0 2 0 0.2 8 0.5 0.75 0 3 0 1 0 >> log.txt
$draw circle_1vortex2.rid 1 1 2 c 2 0 0.15 50 0 -0.25 0 c 3 0 0.2 8 0 0.25 22.5 1 0 -c 0 0.3 >> log.txt

end_time=$(date +%s%3N) 

# convert files to gnuplot:
for f in square1 square2 square3 square4 square5 square6 rectangle1 rectangle2 circle_2vortexes1 circle_2vortexes2 circle_1vortex1 circle_1vortex2
do
	$o2g "$f".rid "$f".ind >> log.txt
	$o2g "$f".rid_nvf_x "$f".ind_nvf_x >> log.txt
	$o2g "$f".rid_nvf_y "$f".ind_nvf_y >> log.txt
done

# gnuplot checking
# paste circle1.ind_nvf_x circle1.ind_nvf_y> circle1.ind_nvf
# plot 'circle1.ind_nvf' every 10:10 u 1:2:($3/50):($7/50)  with vectors
# set pm3d map 
# splot 'circle1' u 1:2:3 


$1 draw.fmm >>log.txt

# An awk script is used to calculate the sum of the absolute values
# of the differences between the results present in the different files.
# Here the meaningful result is on the third column, being the result 
# of a mode calculation (real part).

# as refractive index is 1 or 2 and NbP = 512. 1024 corresponds to 2 line difference.
awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {comp=0; norm=0;} \
	FNR>1&&FNR==NR{a[FNR]=abs(\$3); next} \
	FNR>1&&a[FNR]!=abs(\$3){comp+=abs(abs(\$3)-a[FNR])} \
	END{\
		Nlines=sqrt(FNR);\
		if (comp>4.5*Nlines){ \
			exit 0 \
		}\
		exit 1 ;\
	}"
	
#echo $awk_script > awk_scr.awk 


# tests squares
for f in square1.ind square2.ind square3.ind square4.ind square5.ind square6.ind square7.ind square8.ind square9.ind square10.ind square11.ind square12.ind
do

	if awk -f awk_scr.awk ref_square.ind $f 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m $f \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_square.ind_nvf_x "$f"_nvf_x 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_x \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_square.ind_nvf_y "$f"_nvf_y
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_y \033[0m"
  	fi

done

# tests rectangles
for f in rectangle1.ind rectangle2.ind
do

	if awk -f awk_scr.awk ref_rectangle.ind $f 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m $f \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_rectangle.ind_nvf_x "$f"_nvf_x 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_x \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_rectangle.ind_nvf_y "$f"_nvf_y
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_y \033[0m"
  	fi

done

# tests circle_2vortexes
for f in circle_2vortexes1.ind circle_2vortexes2.ind
do

	if awk -f awk_scr.awk ref_circle_2vortexes.ind $f 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m $f \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_circle_2vortexes.ind_nvf_x "$f"_nvf_x 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_x \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_circle_2vortexes.ind_nvf_y "$f"_nvf_y
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_y \033[0m"
  	fi

done

# tests circle_1vortexes
for f in circle_1vortex1.ind circle_1vortex2.ind
do

	if awk -f awk_scr.awk ref_circle_1vortex.ind $f 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m $f \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_circle_1vortex.ind_nvf_x "$f"_nvf_x 
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_x \033[0m"
  	fi
  	
	if awk -f awk_scr.awk ref_circle_1vortex.ind_nvf_y "$f"_nvf_y
	then
  		test_fail=1
  		printf "\033[1m\x1b[31m "$f"_nvf_y \033[0m"
  	fi

done


elapsed_seconds=$((end_time - start_time))
#echo "duration" $elapsed_seconds

if [ "$test_fail" -eq 1 ] 
then
    printf "\033[1m\x1b[31m FAIL \033[0m\n"
else
  printf " \x1b[32mOK\033[0m\n"
  rm square*
  rm rectangle*
  rm circle*
fi

exit $test_fail
