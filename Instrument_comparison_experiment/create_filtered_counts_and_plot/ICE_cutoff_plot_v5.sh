#This script needs :
#1  paths to folders containing all Rspecies_counts
#2 folder with all the expected counts from mock samples - set name in Mockcounts variable

#also requires the following scripts:
#   R_cutoff_proportional.R script
#   Mock_proportional.R
#   merge_plots_v3.R

# the path to these three directories need to be entered into R_cutoff_proportional.R script too
#set directories of R counts
TEF_MiSeq_Rcounts=""
TEF_iSeq_Rcounts=""
ITS_PE_MiSeq_Rcounts=""
ITS_PE_iSeq_Rcounts=""
ITS_PE_iSeqrep2_Rcounts=""

#record paths for all data used in a text file
echo $TEF_MiSeq_Rcounts > input_data
echo $TEF_iSeq_Rcounts >> input_data
echo $ITS_PE_MiSeq_Rcounts >> input_data
echo $ITS_PE_iSeq_Rcounts >> input_data
echo $ITS_PE_iSeqrep2_Rcounts >> input_data

Mockcounts="Mock_counts_v6"
workdir=$(pwd)
#remove any files from prior analysis
if [ -d "TEF_data" ]
then
  #remove old files ######
  echo "done this before... removing previous analysis..."
  rm -r ITS_data
  rm -r TEF_data

  rm $workdir/"$Mockcounts"/*proportional.csv
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

#copy Rcounts to new folders in working dir
cp $TEF_MiSeq_Rcounts/*Rspecies* TEF_data/TEF_merge/TEF_MiSeq
cp $TEF_iSeq_Rcounts/*Rspecies* TEF_data/TEF_merge/TEF_iSeq
cp $ITS_PE_MiSeq_Rcounts/*Rspecies* ITS_data/ITS_PE/ITS_MiSeq
cp $ITS_PE_iSeq_Rcounts/*Rspecies* ITS_data/ITS_PE/ITS_iSeq
cp $ITS_PE_iSeqrep2_Rcounts/*Rspecies* ITS_data/ITS_PE/ITS_iSeqrep2

#rename count files to include instrument and add instrument to lines (iSeq)
#also take only first two sections of species name
#(ie. genus and first species. If multiple species possible, only first is taken for purpose of matching with expected)
for file in */*/*_iSeq/*Rspecies_counts.csv
do
  path=$( echo $file | cut -d"/" -f1-3)
  newname=$( echo $file | cut -d"/" -f4 | sed 's/.Rspecies/_iSeq_species/' | cut -d"_" -f2-5 )
  mv $file "$path"/"$newname"

  sample=$( echo $newname | cut -d"_" -f1-2)
  cat "$path"/"$newname" | awk -F"," -v where=`basename $sample` '{print $1","$2","where}' >"$path"/"$sample".csv
  rm "$path"/"$newname"
  #cut to first species
  while read p
  do
    cut_species=$(echo $p | cut -d"," -f2 | cut -d"-" -f1)
    echo $p | awk -F"," -v where=`basename $cut_species` '{print $1","where","$3}' >> "$path"/mod_"$sample".csv
  done<"$path"/"$sample".csv
  #remove 'mod_' from filename
  mv "$path"/mod_"$sample".csv "$path"/"$sample".csv
done
#iseq rep2
for file in */*/*_iSeqrep2/*Rspecies_counts.csv
do
  path=$( echo $file | cut -d"/" -f1-3)
  newname=$( echo $file | cut -d"/" -f4 | sed 's/.Rspecies/_iSeqrep2_species/' | cut -d"_" -f2-5 )
  mv $file "$path"/"$newname"

  sample=$( echo $newname | cut -d"_" -f1-2)
  cat "$path"/"$newname" | awk -F"," -v where=`basename $sample` '{print $1","$2","where}' >"$path"/"$sample".csv
  rm "$path"/"$newname"
  #cut to first species
  while read p
  do
    cut_species=$(echo $p | cut -d"," -f2 | cut -d"-" -f1)
    echo $p | awk -F"," -v where=`basename $cut_species` '{print $1","where","$3}' >> "$path"/mod_"$sample".csv
  done<"$path"/"$sample".csv
  #remove 'mod_' from filename
  mv "$path"/mod_"$sample".csv "$path"/"$sample".csv
done

#rename count files to include instrument and add instrument to lines (MiSeq)
for file in */*/*_MiSeq/*Rspecies_counts.csv
do
  path=$( echo $file | cut -d"/" -f1-3)
  newname=$( echo $file | cut -d"/" -f4 | sed 's/.Rspecies/_MiSeq_species/' | cut -d"_" -f2-5 )
  mv $file "$path"/"$newname"

  sample=$( echo $newname | cut -d"_" -f1-2)
  cat "$path"/"$newname" | awk -F"," -v where=`basename $sample` '{print $1","$2","where}' >"$path"/"$sample".csv
  rm "$path"/"$newname"
  #cut to first species
  while read p
  do
    cut_species=$(echo $p | cut -d"," -f2 | cut -d"-" -f1)
    echo $p | awk -F"," -v where=`basename $cut_species` '{print $1","where","$3}' >> "$path"/mod_"$sample".csv
  done<"$path"/"$sample".csv
  #remove 'mod_' from filename
  mv "$path"/mod_"$sample".csv "$path"/"$sample".csv
done

#RUN R SCRIPT - remove species below cutoff threshold and calculate proportions

# /usr/bin/Rscript
cd $workdir/TEF_data/TEF_merge/TEF_MiSeq
Rscript $workdir/R_cutoff_proportional_v2.R
cd $workdir/TEF_data/TEF_merge/TEF_iSeq
Rscript $workdir/R_cutoff_proportional_v2.R
cd $workdir/ITS_data/ITS_PE/ITS_MiSeq
Rscript $workdir/R_cutoff_proportional_v2.R
cd $workdir/ITS_data/ITS_PE/ITS_iSeq
Rscript $workdir/R_cutoff_proportional_v2.R
cd $workdir/ITS_data/ITS_PE/ITS_iSeqrep2
Rscript $workdir/R_cutoff_proportional_v2.R

cd $workdir/$Mockcounts


#generate proportional Mockcounts
Rscript $workdir/Mock_proportional.R


##################### ITS PE merging & plotting #####################
# make new comparisons directory and merge corresponding mock and sample counts files
cd $workdir
mkdir ITS_data/ITS_PE/comparisons
mv ITS_data/ITS_PE/*/*cutoff_filtered_* ITS_data/ITS_PE/comparisons
cp $workdir/"$Mockcounts"/*Mock_proportional.csv ITS_data/ITS_PE/comparisons
cd ITS_data/ITS_PE/comparisons
###
#merge each 4 files (mock,iseq,iseqrep2 and miseq) into one counts file for plotting
mkdir merged

for file in *Mock_proportional.csv
do
pre=$( echo $file | cut -d"." -f1)
cat "$pre"* | grep -v Counts | sed 's/"//g' > merged/"$pre"_merge.csv
done
####

#move into merged folder and plot data #############
cd merged
#RUN R SCRIPT - plot merged counts
# /usr/bin/Rscript
Rscript $workdir/merge_plots_v3.R

##################### TEF merge merging & plotting ######################
cd $workdir
mkdir TEF_data/TEF_merge/comparisons
mv TEF_data/TEF_merge/*/*cutoff_filtered_* TEF_data/TEF_merge/comparisons
cp $workdir/"$Mockcounts"/*Mock_proportional.csv TEF_data/TEF_merge/comparisons
cd TEF_data/TEF_merge/comparisons
###
#merge each 3 files (mock,iseq and miseq) into one counts file for plotting
mkdir merged

for file in *Mock_proportional.csv
do
pre=$( echo $file | cut -d"." -f1)
cat "$pre"* | grep -v Counts | sed 's/"//g' > merged/"$pre"_merge.csv
done

#move into merged folder and plot data ##############
cd merged
#RUN R SCRIPT - plot merged counts
# /usr/bin/Rscript
Rscript $workdir/merge_plots_v3.R

##########################################################
cd $workdir
