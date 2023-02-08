#Instrument comparison experiment analysis code 

#starting point: have run TEF merge and ITS PE pipelines on iSeq and MiSeq data
# Set paths to species count data from mapping pipelines
TEF_MiSeq="<EDIT PATH HERE !>/TEF_merge_analysis/demulti_fastq/Primer_QC/TEF_v12.4_0.06/Counts/Species"
TEF_iSeq="<EDIT PATH HERE !>/TEF_merge_analysis/demulti_fastq/Primer_QC/TEF_v12.4_0.06/Counts/Species"
ITS_PE_MiSeq="<EDIT PATH HERE !>/ITS_PE_analysis/Primer_QC/v2_ISHAM_ITS_0.06/Counts/Species"
ITS_PE_iSeq="<EDIT PATH HERE !>/ITS_PE_analysis/Primer_QC/v2_ISHAM_ITS_0.06/Counts/Species"
ITS_PE_iSeqrep2="<EDIT PATH HERE !>/ITS_PE_analysis/Primer_QC/v2_ISHAM_ITS_0.06/Counts/Species"

## Initial count file prep ####
#Remove any human counts from data, using following code (run in each of the directories containing count files):
#create new count files without human counts -  with prefix "fungal_"
for file in *.csv; do echo $file; grep -v Homo_sapiens $file > fungal_"$file"; done
#create new dir for fungal only data
mkdir fungal_only
#move fungal only count files into new dir
mv fungal* fungal_only/
#enter new dir and remove "fungal_" prefix from filenames
cd fungal_only/
for file in *.csv; do newname=$( echo $file | sed 's/fungal_//' ); mv $file $newname ; done

# Starting in the directory you would like to perform the analysis, and run code below: ####################################################
mkdir R_analysis
cd R_analysis/

## Create cutoff filtered counts and plot data as raw and proportional barplots ####
#       (all data for each sample is merged together) 
./ICE_cutoff_plot_v5.sh
# NOTE! this script and a required script need path variables to be set within them ! ###
# Requires:
#   R_cutoff_proportional.R script ( path variables to be set here too !)
#   Mock_proportional.R
#   merge_plots_v3.R

## Calculate Identification rates ##########################################################################################################
cd ..
mkdir ID_rate
cd ID_rate

#calculate the ID rates and output a matrix
./species_IDrates_v4.sh
#   Requires: R_count_matrix_v3.R

# use matrix of ID rates to plot heatmaps 
Rscript plot_IDrate_heatmap_v3.R

## Calculate ID rates for human background experiment ##
cd ..
mkdir TEF_human_background_comparison
cd TEF_human_background_comparison
mkdir TEF_no_human
mkdir TEF_human_10X
cd TEF_no_human
#calculate the ID rates and output a matrix
./species_IDrates_v4.sh

cd ../TEF_human_10X/
cp ../TEF_no_human/species_IDrates_v4.sh .
#calculate the ID rates and output a matrix
./species_IDrates_v4.sh

cd ..
#plot ID rates for human background experiment (TEF) 
Rscript R_plotting_human_v_nonhuman_TEF_v3.R

cd ..
mkdir ITS_human_background_comparison
cd ITS_human_background_comparison/
mkdir ITS_no_human
mkdir ITS_human_10X
cd ITS_no_human/
#calculate the ID rates and output a matrix
./species_IDrates_v4.sh

cp species_IDrates_v4.sh ../ITS_human_10X/
cp R_count_matrix_v3.R ../ITS_human_10X/

cd ../ITS_human_10X/
#calculate the ID rates and output a matrix
./species_IDrates_v4.sh
cd ..
#plot ID rates for human background experiment (ITS1)
Rscript R_plotting_human_v_nonhuman_ITS1_v3.R

cd ..

### Compare expected species with experimental results ###########################################################################
mkdir Expected_species_only
cd Expected_species_only

./expectedspecies_analysis_v6.sh
#Requires:
#   adjust_proportions.R

./calc_chisq_all_samples_v4.sh


cd ..
mkdir species_quant
cd species_quant
#Add target name to sample column and create merged csv files for all samples:
./Species_quant_v2.sh
#### The resulting merged files can then be used for any further analyses of subsets of samples #############
#     Eg. for performing comparison of species quantitation by each instrument, 
#       gather relevant samples and use the following:  
#first check if data is normal or not 
TEF_v_ITS_normality_kruskal.R
#perform comparison 
Instrument_comparison_stats_v3.R

# For human bakground comparison instead, use the following script: 
Human_bg_comparison_stats_v2.R

#### Dendrogram & barplotting for example mock vs. experimental community figures ################################################
mkdir dendro_barplotting
cd dendro_barplotting
##MA1
cp $TEF_MiSeq/*MA1.[1-2].unique_species .

for file in TEF*[0-9]_MA*unique_species
do
  new=$( echo $file | awk -F "_" '{print$1"_MiSeq_"$2"_"$3}')
  mv $file $new
done


cp $TEF_iSeq/*MA1.[1-2].unique_species .

for file in TEF*[0-9]_MA*unique_species
do
  new=$( echo $file | awk -F "_" '{print$1"_iSeq_"$2"_"$3}')
  mv $file $new
done


cp $ITS_PE_MiSeq/*MA1.[1-2]_unique_species .
for file in ITS*[0-9]_MA*unique_species
do
  new=$( echo $file |  awk -F "_" '{print$1"_MiSeq_"$2"_"$3"_"$4}')
  mv $file $new
done


cp $ITS_PE_iSeq/*MA1.[1-2]_unique_species .

for file in ITS*[0-9]_MA*unique_species
do
  new=$( echo $file |  awk -F "_" '{print$1"_iSeq_"$2"_"$3"_"$4}')
  mv $file $new
done

cp $ITS_PE_iSeqrep2/*MA1.[1-2]_unique_species .

for file in ITS*[0-9]_MA*unique_species
do
  new=$( echo $file |  awk -F "_" '{print$1"_iSeqrep2_"$2"_"$3"_"$4}')
  mv $file $new
done


for file in *_unique_species; do new=$( echo $file | sed 's/_unique_species/.unique_species/' ); mv $file $new; done

./join_TEF_v2.sh

#Create the plots 
#NOTE - need to change the "samplename" variable to 'MA1'
Rscript phyloseq_dendro_barplot_v3.R
#    Requires: fungi_tax_TEF_v4.tsv

##MD1

cd ../MD1

cp $TEF_MiSeq/*MD1.[1-2].unique_species .

for file in TEF*[0-9]_MD*unique_species
do
  new=$( echo $file | awk -F "_" '{print$1"_MiSeq_"$2"_"$3}')
  mv $file $new
done


cp $TEF_iSeq/*MD1.[1-2].unique_species .

for file in TEF*[0-9]_MD*unique_species
do
  new=$( echo $file | awk -F "_" '{print$1"_iSeq_"$2"_"$3}')
  mv $file $new
done


cp $ITS_PE_MiSeq/*MD1.[1-2]_unique_species .
for file in ITS*[0-9]_MD*unique_species
do
  new=$( echo $file |  awk -F "_" '{print$1"_MiSeq_"$2"_"$3"_"$4}')
  mv $file $new
done


cp $ITS_PE_iSeq/*MD1.[1-2]_unique_species .

for file in ITS*[0-9]_MD*unique_species
do
  new=$( echo $file |  awk -F "_" '{print$1"_iSeq_"$2"_"$3"_"$4}')
  mv $file $new
done

cp $ITS_PE_iSeqrep2/*MD1.[1-2]_unique_species .

for file in ITS*[0-9]_MD*unique_species
do
  new=$( echo $file |  awk -F "_" '{print$1"_iSeqrep2_"$2"_"$3"_"$4}')
  mv $file $new
done


for file in *_unique_species; do new=$( echo $file | sed 's/_unique_species/.unique_species/' ); mv $file $new; done

./join_TEF_v2.sh

#Create the plots 
#NOTE - need to change the "samplename" variable to 'MD1'
Rscript phyloseq_dendro_barplot_v3.R
#    Requires: fungi_tax_TEF_v4.tsv

#### Plotting all data for each single species samples ############################################################################

#Starting with TEF and ITS filtered data, need to ensure that target name is in filename and samplename
#   Can add this using following code: 
## For TEF:
#Miseq
for file in *_MiSeq_cutoff_filtered_0_2.csv
do
 new=$( echo $file | sed 's/_MiSeq/_TEF_MiSeq/' )
 cat $file | sed 's/_MiSeq/_TEF_MiSeq/' > $new
done
#iSeq
for file in *_iSeq_cutoff_filtered_0_2.csv
do
 new=$( echo $file | sed 's/_iSeq/_TEF_iSeq/' )
 cat $file | sed 's/_iSeq/_TEF_iSeq/' > $new
done

#move into new dir
mkdir TEF_cutoffs/
mv *_TEF_*.csv TEF_cutoffs/

## and for ITS:
#Miseq
for file in *_MiSeq_cutoff_filtered_0_2.csv
do
 new=$( echo $file | sed 's/_MiSeq/_ITS1_MiSeq/' )
 cat $file | sed 's/_MiSeq/_ITS1_MiSeq/' > $new
done

#iSeq
for file in *_iSeq_cutoff_filtered_0_2.csv
do
 new=$( echo $file | sed 's/_iSeq/_ITS1_iSeq/' )
 cat $file | sed 's/_iSeq/_ITS1_iSeq/' > $new
done

mkdir ITS_cutoffs/
mv *_ITS1_*.csv ITS_cutoffs/

###copy TEF, ITS and mocks species samples data to new dir for merging (MANUALLY)

#In the dir containing all the data to be merged, make merged counts for plotting:
cd /home/mfbx9dw5/Dropbox/Instrument_comparison_exp_data/Single_species_samples/R_cutoff_method/combi_single_species
mkdir merged
#merge by sample type
for file in *Mock_proportional.csv; do pre=$( echo $file | cut -d"." -f1); cat "$pre"* | grep -v Counts | sed 's/"//g' > merged/"$pre"_merge.csv; done
cd merged/

#remove mock lines from merged csv files
for file in *_merge.csv; do cat $file | grep -v Mock > mod_"$file"; done
mkdir no_mock
mv mod_* no_mock/
cd no_mock/
for file in mod_*; do new=$( echo $file | sed 's/mod_//'); mv $file $new; done

## Plot single species  using:
barplot_combi_SS_v1.R

#### Plotting limit of detection (LOD) experiment data ############################################################################
# Begin in dir containing relevant LOD sample files (cutoff counts)

#remove mock lines and "Miseq" or "iSeq" from merged csv files
for file in *_merge.csv; do cat $file | grep -v Mock | grep -v _MiSeq  | grep -v _iSeq > mod_"$file"; done
mkdir no_mock
mv mod_* no_mock/
cd no_mock/
for file in mod_*; do new=$( echo $file | sed 's/mod_//'); mv $file $new; done

mkdir merged
#merge by sample type
for file in *Mock_proportional.csv; do pre=$( echo $file | cut -d"." -f1); cat "$pre"* | grep -v Counts | sed 's/"//g' > merged/"$pre"_merge.csv; done
cd merged/

#remove mock lines from merged csv files
for file in *_merge.csv; do cat $file | grep -v Mock > mod_"$file"; done
for file in mod_*; do new=$( echo $file | sed 's/mod_//'); mv $file $new; done

#plot data using: 
barplot_combi_LOD_v1.R
