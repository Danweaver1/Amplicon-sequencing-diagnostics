#!/bin/bash
#$ -cwd
#$ -V
#joins .unique_species files into a matrix
if [ -e out ]
then
  rm out
  rm *.join
  rm sorted_*
fi


for file in *unique_species
do
        if [ -e $file ]
        then
		PREFIX=`echo $file | sed -e 's/\.unique_species//'`
		OUTFILE="$PREFIX.join"
    cat $file | awk '{$1=$1};1' | sort -t" " -k2 > sorted_"$PREFIX"
		join --nocheck-order -1 1 -2 2 -a 1 -e '0' -o '0,2.1' TEFv4_list2 sorted_"$PREFIX" > $OUTFILE
		echo -e ' '"$OUTFILE"'\n'"$(cat -- "$OUTFILE")" > "$OUTFILE"
		sed -i 's/\.join//g' $OUTFILE
		paste TEFv4_list3 *.join  | sed 's/\s/\t/g' |\
		cut -f1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,70 > out

	else
                echo "join couldn't run, please check you have the right files"
        fi
done
