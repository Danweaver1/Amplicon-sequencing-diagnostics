
## Calculate identification rates for each species present at a low level in the mock community ###
#     (ie. dominant taxa for each mock is ignored)
#Requirements:
# R_count_matrix_v3.R
# must provide files with list of samples to analyse for each target named "TEF_sample_set" and "ITS_sample_set"
# Set the directory for proportional mock counts ( not including replicates)
mockcounts_dir="<PATH>/Mock_counts_v6"
R_analysis_dir=$( dirname $mockcounts_dir )
# R cutoff counts should be within same directory as mock counts, with dir layout as below:
TEF_Rcounts=""$R_analysis_dir"/TEF_data/TEF_merge/comparisons"
ITS_PE_Rcounts=""$R_analysis_dir"/ITS_data/ITS_PE/comparisons"

Mockcounts=$( basename $mockcounts_dir)
if [ -d "$Mockcounts" ]
then
  rm "$Mockcounts"/*
else
  mkdir "$Mockcounts"
fi

cp $mockcounts_dir/*Mock_proportional* $Mockcounts

workdir=$(pwd)
#remove any files from prior analysis
if [ -d "TEF_data" ]
then
  rm -r ITS_data
  rm -r TEF_data
  rm -r headerless_counts
  rm All_R_counts.csv
  rm ID_species_R_counts.csv
  rm Missing_species_R_counts.csv
  rm ID_rate.csv
fi

#create organised folders for each amplicon and sequencer combination
mkdir ITS_data
mkdir ITS_data/ITS_PE
mkdir ITS_data/ITS_PE/ITS_MiSeq
mkdir ITS_data/ITS_PE/ITS_iSeq
mkdir ITS_data/ITS_PE/ITS_iSeqrep2

mkdir TEF_data
mkdir TEF_data/TEF_merge
mkdir TEF_data/TEF_merge/TEF_iSeq
mkdir TEF_data/TEF_merge/TEF_MiSeq

#copy Rcounts to new folders in working dir
cp $TEF_Rcounts/*_MiSeq_cutoff_filtered_* TEF_data/TEF_merge/TEF_MiSeq
cp $TEF_Rcounts/*_iSeq_cutoff_filtered_* TEF_data/TEF_merge/TEF_iSeq
cp $ITS_PE_Rcounts/*_MiSeq_cutoff_filtered_* ITS_data/ITS_PE/ITS_MiSeq
cp $ITS_PE_Rcounts/*_iSeq_cutoff_filtered_* ITS_data/ITS_PE/ITS_iSeq
cp $ITS_PE_Rcounts/*_iSeqrep2_cutoff_filtered_* ITS_data/ITS_PE/ITS_iSeqrep2

#cp  $TEF_iSeq_Rcounts TEF_data/TEF_merge

#make headerless mock counts files necessary for analysis using sample set files
# also make headerless mock counts with lines for 1or2% species only
for file in *sample_set
  do
    while read p
    do
      #take sample name (don't include rep number) from sample set file line
      sample=$( echo $p | cut -d"." -f1 )
      #find mock counts file
      for f in $Mockcounts/"$sample".[0-9]_Mock_proportional*.csv
      do
        fnopath=$( echo $f | cut -d"/" -f2 )
        #make counts files without header (in workdir)
        #cat $f | grep -v Counts > headerless_"$fnopath"
        #for analysing low species only, add # to beginning of previous line and remove # from next line
        cat $f | grep -v Counts |  grep -E '(0.02|0.01)' > headerless_"$fnopath"

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
  for f in headerless_"$sample".[0-9]_Mock_proportional*.csv
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
        #if species found/not found in counts file, add to ID_species/Missing_species file in corresponding dir
        if grep -Fq "$species" "$f2"
        then
          echo "$species" >> $filepath/ID_species
        else
          echo "$species" >> $filepath/Missing_species
        fi
      done
    done<$f
  done

done<TEF_sample_set

#analyse ITS
while read p
do
  #using sample set to dictate which mock counts to look into
  #take sample name (don't include rep number) from sample set file line
  sample=$( echo $p | cut -d"." -f1 )
  echo $sample
  #find mock counts file
  for f in headerless_"$sample".[0-9]_Mock_proportional*.csv
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
        #if species found/not found in counts file, add to ID_species/Missing_species file in corresponding dir
        if grep -Fq "$species" "$f2"
        then
          echo "$species" >> $filepath/ID_species
        else
          echo "$species" >> $filepath/Missing_species
        fi
      done
    done<$f
  done

done<ITS_sample_set

mkdir headerless_counts
mv headerless_*.csv headerless_counts


#generate ID/missing counts files
for file in */*/*/ID_species
do
  filepath=$( echo $file | cut -d"/" -f1-3 )
  #generate standard counts in each dir
  cat $file | sort | uniq -c | sort -nr > $filepath/ID_species_counts &&
  #generate global counts in format for R with pipeline and instrument info
  PREFIX=$( echo $filepath | cut -d"/" -f2-3 | sed -e 's/\//_/' | cut -d"_" -f2-4 | awk '{ print $0"_ID"}')
  cat $filepath/ID_species_counts | sed -e 's/^ *//;s/ /,/' | awk -v where=`basename $PREFIX` '{print $1","where}' >> ID_species_R_counts.csv

done

#generate ID/missing counts files
for file in */*/*/Missing_species
do
  filepath=$( echo $file | cut -d"/" -f1-3 )
  #generate standard counts in each dir
  cat $file | sort | uniq -c | sort -nr > $filepath/Missing_species_counts &&
  #generate global counts in format for R with pipeline and instrument info
  PREFIX=$( echo $filepath | cut -d"/" -f2-3 | sed -e 's/\//_/' | cut -d"_" -f2-4 | awk '{ print $0"_Missing"}')
  cat $filepath/Missing_species_counts | sed -e 's/^ *//;s/ /,/' | awk -v where=`basename $PREFIX` '{print $1","where}' >> Missing_species_R_counts.csv

done

#if Missing file is empty, fill contents with empty equivalent of found
if [ -s "Missing_species_R_counts.csv" ]
then
  cat *species_R_counts.csv > All_R_counts.csv
else
  cat ID_species_R_counts.csv | sed 's/[0-9]*,/0,/' | sed 's/_ID/_Missing/' > Missing_species_R_counts.csv &&
  cat *species_R_counts.csv > All_R_counts.csv

fi
#then run R script to convert the ID rates into a matrix 
Rscript R_count_matrix_v3.R
