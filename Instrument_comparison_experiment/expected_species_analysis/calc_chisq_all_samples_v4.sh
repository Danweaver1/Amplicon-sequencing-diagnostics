# Run chi squared tests on data from instrument comparison experiment ####

#NB mock counts need to include seperate files for replicate samples for this analysis to assess them!!
#Requires:
#   chi_iSeq_v3.R
#   chi_MiSeq_v2.R
#Need to add paths to proportional counts files: 
TEF_iSeq_data=""
ITS_iSeq_data=""
TEF_MiSeq_data=""
ITS_MiSeq_data=""
Mock_data=""

#remove any files from prior analysis
if [ -d "iSeq_chi" ]
then
  rm -r iSeq_chi/
  rm -r MiSeq_chi/
fi


### iSeq ##
mkdir iSeq_chi
cp $TEF_iSeq_data/TEF_Exp_sp*prop_adjusted.csv iSeq_chi
cp $ITS_iSeq_data/ITS_Exp_sp*prop_adjusted.csv iSeq_chi
#copy replicate mock count data
cp $Mock_data/M*proportional.csv iSeq_chi

cd iSeq_chi
#Perform chi squared tests
Rscript ../chi_iSeq_v3.R
mkdir p_vals
mv *vals.csv p_vals
cd ..

### Miseq ##
mkdir MiSeq_chi
cp $TEF_MiSeq_data/TEF_Exp_sp*prop_adjusted.csv MiSeq_chi
cp $ITS_MiSeq_data/ITS_Exp_sp*prop_adjusted.csv MiSeq_chi
cp $Mock_data/M*proportional.csv MiSeq_chi

cd MiSeq_chi

#Perform chi squared tests
Rscript ../chi_MiSeq_v2.R
mkdir p_vals
mv *vals.csv p_vals
