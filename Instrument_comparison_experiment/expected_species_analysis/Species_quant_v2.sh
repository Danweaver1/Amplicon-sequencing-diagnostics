
#remove any files from prior analysis
if [ -d "unmod_counts" ]
then
  rm -r unmod_counts/
  rm -r merged/
  rm -r mod_counts/
  rm -r mock_counts/
fi

#gather expected species only counts with adjusted proportions
##sample data
cp <EDIT PATH HERE !>/Expected_species_only/MiSeq_chi/*prop_adjusted.csv .
cp <EDIT PATH HERE !>/Expected_species_only/iSeq_chi/*prop_adjusted.csv .
#mock data - includes reps
mkdir mock_counts
cp <EDIT PATH HERE !>/Expected_species_only/MiSeq_chi/*Mock_proportional.csv .


#add target info to Sample field
#ITS
for file in ITS*iSeq_*
do
  cat $file | sed 's/iSeq/iSeq_ITS/' > mod_"$file"
done

#ITS iseqrep2
for file in ITS*iSeqrep2_*
do
  cat $file | sed 's/iSeqrep2/iSeqrep2_ITS/' > mod_"$file"
done

#TEF
for file in TEF*iSeq_*
do
  cat $file | sed 's/iSeq/iSeq_TEF/' > mod_"$file"
done


mkdir unmod_counts
mv ITS* unmod_counts
mv TEF* unmod_counts

mkdir merged

for file in *Mock_proportional.csv
do
pre=$( echo $file | cut -d"_" -f1)
cat *"$pre"* | grep -v Counts | sed 's/"//g' > merged/"$pre"_merge.csv
done

mkdir mod_counts
mv mod_*csv mod_counts
mv *Mock_proportional.csv mock_counts/
