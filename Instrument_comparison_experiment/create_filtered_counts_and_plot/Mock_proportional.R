### Generate proportional counts of Mock count data ####
###This script takes in expected counts with following columns:
#   Counts,Species,Sample
# and adds a new column containing counts converted into proportions

#### generate proportional counts in R ####
library("stringr")
full_file_list=(list.files(pattern=".csv"))
file_list=full_file_list[file.size(full_file_list) > 0]
library("dplyr")
for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,".csv","")
  #read counts into R
  counts=as.data.frame(read.csv(filename, header=TRUE, sep=",", col.names=c("Counts","Species","Sample")))
  #add a Proportional_counts column
  counts[4]=counts["Counts"]/sum(counts["Counts"])
  colnames(counts)[4] = "Proportional_counts"
  write.csv(counts, paste(samplename,"_proportional.csv", sep=""), row.names=F) #create a file of full counts
}
