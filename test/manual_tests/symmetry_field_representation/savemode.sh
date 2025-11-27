# save all picture corresponding to the first modes
# syntaxe: savemode.sh $1 $2
# first input parameter $1 = numberofmodes
# second input parameter $2 = filename

NBMODES=$1
FILE=$2
#save pictures:
for (( i=1 ; i <= $NBMODES ; i+=1))
do
	# we have to take the ith value (the ith mode)
	j=`expr $i - 1`

	if test $j -le 0
	then
		cut -d\  -f 7 $FILE.out | sed '2,50d'> effectiveindex.txt
		cut -d\  -f 4 $FILE.out | sed '2,50d'| cut -d , -f 1 > modenumber.txt
	else
		string='1,'$j'd'
		cut -d\  -f 7 $FILE.out | sed $string | sed '2,50d' > effectiveindex.txt
		cut -d\  -f 4 $FILE.out | sed $string | sed '2,50d'| cut -d , -f 1 > modenumber.txt
	fi

	# read the two files:
	read effectiveindex < effectiveindex.txt
	read modenumber < modenumber.txt

	echo "$i saving mode $modenumber with $effectiveindex"

	echo "	set pm3d map
			set terminal png size 800,800
			set output 'Ex_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Ex_m_$modenumber.mode' 
			set output 'Ey_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Ey_m_$modenumber.mode'
			set output 'Ez_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Ez_m_$modenumber.mode'

			set output 'Ex2_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Ex2_m_$modenumber.mode' 
			set output 'Ey2_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Ey2_m_$modenumber.mode'
			set output 'Ez2_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Ez2_m_$modenumber.mode'

			set output 'Dx_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Dx_m_$modenumber.mode' 
			set output 'Dy_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Dy_m_$modenumber.mode'
			set output 'Dz_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Dz_m_$modenumber.mode'

			set output 'Hx_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Hx_m_$modenumber.mode' 
			set output 'Hy_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Hy_m_$modenumber.mode'
			set output 'Hz_m_$effectiveindex-$FILE-$modenumber.mode.png'
			splot '$FILE._Hz_m_$modenumber.mode'" > plot_mode.plt
	gnuplot plot_mode.plt
done

#clean
rm plot_mode.plt
rm effectiveindex.txt
rm modenumber.txt
