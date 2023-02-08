## Adjusts proportions of species after non-expected species (those not in mock) have been filtered out of the data
library("stringr")
full_file_list=(list.files(pattern='.*Exp_sp'))
file_list=full_file_list[file.size(full_file_list) > 0]
library("dplyr")
for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,".csv","")
  #read counts into R
  counts=as.data.frame(read.csv(filename, header=TRUE, sep=","))
  #adjust Proportional_counts column
  if (sum(counts[1] > 0)) {
    counts[4]=counts["Counts"]/sum(counts["Counts"])
  }
  write.csv(counts, paste(samplename,"_prop_adjusted.csv", sep=""), row.names=F) #create a file of full counts
}

