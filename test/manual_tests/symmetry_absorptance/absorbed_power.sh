# extract the absorbed power from afmm output
# syntaxe: absorbed_power.sh $1 $2
# first input parameter $1 = afmm_output_filename
# second input parameter $2 = report_filename

FILE=$1
REPORT=$2
#extract absorbed power


# extract the number of colomn in the file
#awk '{print NF}' $FILE.out > nbcol.txt
#echo "nbcol : "
#cat nbcol.txt
# for each row, extract the last colomn:
echo "" > absorptance.txt
while read line  
do 
	# extract the number of colomn 
	echo $line | awk '{print NF}' > nbcol.txt
	read nbcol < nbcol.txt

	# extract the value of the last colomn
    echo $line | cut -d\  -f $nbcol >> absorptance.txt

done < $FILE.out

# extract the absorbed power
# generation
sed '2,9d' absorptance.txt > absorptance2.txt
read generation < absorptance2.txt
# generation averaged over a radius
sed '1,2d' absorptance.txt | sed '2,8d' > absorptance2.txt
read generation_averaged < absorptance2.txt
# power
sed '1,3d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power1 < absorptance2.txt
sed '1,4d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power2 < absorptance2.txt
sed '1,5d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power3 < absorptance2.txt
python -c "print $power1-$power2-$power3" > absorptance2.txt
read power < absorptance2.txt



# powerZ
sed '1,6d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power1 < absorptance2.txt
sed '1,7d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power2 < absorptance2.txt
python -c "print $power1-$power2" > absorptance2.txt
read powerZ < absorptance2.txt
# monitor
sed '1,8d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power1 < absorptance2.txt
sed '1,9d' absorptance.txt | sed '2,8d' > absorptance2.txt
read power2 < absorptance2.txt
python -c "print $power1-$power2" > absorptance2.txt
read monitor < absorptance2.txt

# echo "filename generation generation power powerZ monitor"
echo $FILE $generation $generation_averaged $power $powerZ $monitor
echo $FILE $generation $generation_averaged $power $powerZ $monitor >> $REPORT

# clean
rm absorptance.txt
rm absorptance2.txt
rm nbcol.txt
