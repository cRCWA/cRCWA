#!/bin/bash

test_fail=0

if [[ -z "$1" ]]
then
	printf "$0 should be called referring to AFMM executable.\n"
	printf "For example: $0 ../../afmm\n"
	exit 1
fi
printf "Testing monitor 2:           "

$1 monitor.fmm >log.txt


# An awk script is used to compare the power computed by test and reference
awk_script="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {error=0;} \
	FNR==NR{id[FNR]=\$1; p[FNR]=\$2; pz[FNR]=\$3; m1[FNR]=\$4; m2[FNR]=\$5; m3[FNR]=\$6; next} \
	abs(p[\$1]-\$2)>1e-18 || abs(pz[\$1]-\$3)>1e-18 || abs(m1[\$1]-\$4)>1e-18 || abs(m2[\$1]-\$5)>1e-18 || abs(m3[\$1]-\$6)>1e-18 \
	{error=1; print \$1\" \" p[\$1] \" \" pz[\$1] \" \" m1[\$1]\" \"m2[\$1] \" \" m3[\$1] \" vs \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$6} \
	END{exit error }"
	
echo $awk_script > awk_scr.awk 

if awk -f awk_scr.awk test_output.txt test_output_reference.txt 
then
  power_fail=0
  rm Ex.txt
  rm Ey.txt
  rm Hx.txt
  rm Hy.txt
else
  power_fail=1
    printf "\033[1m\x1b[31m power \033[0m\n"
fi

# concatenate files with computed fields:
paste cas6_Ex.txt cas6_Hy.txt cas6_Ey.txt cas6_Hx.txt cas6_Pz.txt > cas6_Ex_Hy_Ey_Hx_Pz.txt
# Pz is computed from Ex*Hy+Ey*Hx and compared with Pz computed from outdata sz:
# Visual check can be obtained with gnuplot:
# plot 'cas6_Ex_Hy_Ey_Hx_Pz.txt' u 1:($4*$9-$5*$10-$14*$19+$15*$20) w p, 'cas6_Ex_Hy_Ey_Hx_Pz.txt' u 1:24 w l 

awk_script_fields="function abs(x){ \
		return ((x < 0.0) ? -x : x)\
	} \
	BEGIN {comp=0; norm=0;} \
	FNR==NR{sz=(\$4*\$9-\$5*\$10-\$14*\$19+\$15*\$20); Pz=\$24 ;norm+=\$24*\$24; comp+=((sz-Pz)*(sz-Pz)) ; next} \
	END{\
		norm=sqrt(norm/FNR);\
		comp=sqrt(comp/FNR);\
		if (comp>1e-3*norm){ \
			exit 1 \
		}\
	}"
	
echo $awk_script_fields > awk_scr_field.awk 


if awk -f awk_scr_field.awk cas6_Ex_Hy_Ey_Hx_Pz.txt
then
  field_fail=0 
  rm cas6_*.txt
else
  field_fail=1
  printf "\033[1m\x1b[31m field \033[0m\n"
fi

if [ "$power_fail" -ne 1 ] && [ "$field_fail" -ne 1 ]
then
	printf " \x1b[32mOK\033[0m\n"	
else
	test_fail=1
	printf "\033[1m\x1b[31m FAIL \033[0m\n"
fi

exit $test_fail
