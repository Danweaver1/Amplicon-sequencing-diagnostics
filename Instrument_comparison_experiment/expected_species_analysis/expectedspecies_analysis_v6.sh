
#must provide files with list of samples to analyse named "sample_set"
#provide directories to find cutoff R counts for ITS & TEF
#provide directory for proportional mock counts (not including reps)
#Requires:
#   adjust_proportions.R

workdir=$(pwd)
#remove any files from prior analysis
if [ -d "TEF_data" ]
then
  rm -r ITS_data
  rm -r TEF_data
fi

mkdir ITS_data
mkdir ITS_data/ITS_PE
mkdir ITS_data/ITS_PE/ITS_MiSeq
mkdir ITS_data/ITS_PE/ITS_iSeq
mkdir ITS_data/ITS_PE/ITS_iSeqrep2

mkdir TEF_data
mkdir TEF_data/TEF_merge
mkdir TEF_data/TEF_merge/TEF_iSeq
mkdir TEF_data/TEF_merge/TEF_MiSeq

#set directories of R cutoff counts
TEF_Rcounts="" # SET DIR HERE ! #####
ITS_PE_Rcounts="" # SET DIR HERE ! #####
#copy Rcounts to new folders in working dir
cp $TEF_Rcounts/*_MiSeq_cutoff_filtered_* TEF_data/TEF_merge/TEF_MiSeq
cp $TEF_Rcounts/*_iSeq_cutoff_filtered_* TEF_data/TEF_merge/TEF_iSeq
cp $ITS_PE_Rcounts/*_MiSeq_cutoff_filtered_* ITS_data/ITS_PE/ITS_MiSeq
cp $ITS_PE_Rcounts/*_iSeq_cutoff_filtered_* ITS_data/ITS_PE/ITS_iSeq
cp $ITS_PE_Rcounts/*_iSeqrep2_cutoff_filtered_* ITS_data/ITS_PE/ITS_iSeqrep2
##create mock counts dir
Mockcounts="Mock_counts_v6"

if [ -d "$Mockcounts" ]
then
  echo "$Mockcounts already exists, removing previous files"
  rm $Mockcounts/*
else
  mkdir "$Mockcounts"
fi


#get approrpiate mock count files based on sample input files

mockcounts_dir=""  # SET DIR HERE ! #####

for file in sample_set
do
  while read p
  do
    #take sample name (don't include rep number) from sample set file line
    sample=$( echo $p | cut -d"." -f1 )
    #find mock counts file
    for f in $mockcounts_dir/"$sample".[0-9]_Mock_proportional.csv
    do
      #copy mock count files to work with
      cp $f $workdir
    done
  done<$file
done


#analyse TEF
while read p
do
  #using sample set to dictate which mock counts to look into
  #take sample name (don't include rep number) from sample set file line
  sample=$( echo $p | cut -d"." -f1 )
  echo $sample
  #find mock counts file
  for f in "$sample".[0-9]_Mock_proportional.csv
  do
    #read headerless mock counts file lines to find species to look for
    while read q
    do
      #take species name from mock count line
      species=$( echo $q | cut -d"," -f2 |sed 's/"//g')
      echo $species

      #search for species name in all sample count file(s)
      for f2 in TEF_data/TEF_*/TEF_*/"$p"*_cutoff_filtered_0_2.csv
      do
        #get file path to seperate results from each dir
        filepath=$( echo $f2 | cut -d"/" -f1-3 )
        filename=$( echo $f2 | cut -d"/" -f4 )
        sample_instrument=$( echo $filename | cut -d"_" -f1-2 )
        #if species found in counts file, add the line hit to new expected only counts file in corresponding dir. If not found, add a line with 0 counts
        if grep -Fq "$species" "$f2"
        then
          grep -F "$species" "$f2" >> $filepath/TEF_Exp_sp_"$filename"
        else
            echo "0,"$species","$sample_instrument",0" >> $filepath/TEF_Exp_sp_"$filename"
        fi
      done
    done<$f
  done

done<sample_set

echo "TEF analysed"
#analyse ITS
while read p
do
  #using sample set to dictate which mock counts to look into
  #take sample name (don't include rep number) from sample set file line
  sample=$( echo $p | cut -d"." -f1 )
  echo $sample
  #find mock counts file
  for f in "$sample".[0-9]_Mock_proportional.csv
  do
    #read headerless mock counts file lines to find species to look for
    while read q
    do
      #take species name from mock count line
      species=$( echo $q | cut -d"," -f2 |sed 's/"//g')
      echo $species

      #search for species name in all sample count file(s)
      for f2 in ITS_data/ITS_*/ITS_*/"$p"*_cutoff_filtered_0_2.csv
      do
        #get file path to seperate results from each dir
        filepath=$( echo $f2 | cut -d"/" -f1-3 )
        filename=$( echo $f2 | cut -d"/" -f4 )
        sample_instrument=$( echo $filename | cut -d"_" -f1-2 )
        #if species found in counts file, add the line hit to new expected only counts file in corresponding dir. If not found, add a line with 0 counts
        if grep -Fq "$species" "$f2"
        then
          grep -F "$species" "$f2" >> $filepath/ITS_Exp_sp_"$filename"
        else
          echo "0,"$species","$sample_instrument",0" >> $filepath/ITS_Exp_sp_"$filename"
        fi
      done
    done<$f
  done

done<sample_set


echo "ITS analysed"
echo "adjusting proportions on expected species count data"
cd ITS_data/ITS_PE/ITS_MiSeq
Rscript ../../../adjust_proportions.R
cd ../ITS_iSeq
Rscript ../../../adjust_proportions.R
cd ../ITS_iSeqrep2
Rscript ../../../adjust_proportions.R
cd ../../../TEF_data/TEF_merge/TEF_iSeq
Rscript ../../../adjust_proportions.R
cd ../TEF_MiSeq
Rscript ../../../adjust_proportions.R

cd ../../../
mv *_Mock_proportional*.csv $Mockcounts
