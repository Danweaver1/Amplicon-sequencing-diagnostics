###This script takes in expected and sample counts with filenames as "<samplename>_Mock/iSeq/MiSeq.csv" (eg. S05.1_Mock.csv) with format(no header):
#Counts,Species,Sample
  
#### generate cutoff proportional counts in R ####    
# removes species below a proportioncutoff (decimal) set in 'val'
#   eg. val =0.002 removes any species present at 2% or lower
library("stringr")
full_file_list=(list.files(pattern=".csv"))
file_list=full_file_list[file.size(full_file_list) > 0]
library("dplyr")
for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,".csv","")
  #read counts into R
  counts=as.data.frame(read.csv(filename, header=FALSE, sep=",", col.names=c("Counts","Species","Sample")))
  #add a Proportional_counts column
  counts[4]=counts["Counts"]/sum(counts["Counts"])
  colnames(counts)[4] = "Proportional_counts"
  val=0.002
  #calculate cutoff as a percentage (for labelling filenames)
  percent_cutoff=val*100
  percent_cutoff_nopunc=sub("[.]","_",percent_cutoff)
  #filter counts by cutoff
  cutoff_filter_counts=counts %>% filter(Proportional_counts > val)
  cutoff_filter_counts[4]=cutoff_filter_counts["Counts"]/sum(cutoff_filter_counts["Counts"])
  #create a file of cutoff filtered counts - include percentage cutoff used in the filename
  write.csv(cutoff_filter_counts, paste(samplename,"_cutoff_filtered_",percent_cutoff_nopunc,".csv", sep=""), row.names=F)
  counts[5]=counts[4]
  colnames(counts)[5] = "Proportional_counts_cutoff_filtered"
  counts[counts$Proportional_counts_cutoff_filtered < val, "Proportional_counts_cutoff_filtered"] = NA
  write.csv(counts, paste(samplename,"_proportional.csv", sep=""), row.names=F) #create a file of full counts
}






#########
  ####### merge relevant count data with expected counts in bash using merge_cutoff_proportional_counts.sh ###########

#for file in 
#do
#take sample name (not including rep info)
#pre=$( echo $file | cut -d"_" -f2 | cut -d"." -f1)
#cat ITS*"$pre"* | grep -v Counts | sed 's/"//g' > "$pre"_merge.csv
#done

    

#
#
#
#
#


