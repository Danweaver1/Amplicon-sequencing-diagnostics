######input identified and missing species counts

ID_species = "All_R_counts.csv"

MC=read.csv(ID_species, header=FALSE, sep=",",col.names=c("Counts","Species","Sample"))

library(tidyr)
data_wide <- spread(MC, Sample, Counts)
#make a dataframe
df <- as.data.frame(data_wide)
#change all NA to zero
df[is.na(df)] <- 0
#use species names as rownames
rownames(df) <- df[,1]
#remove species column
df[,1] <- NULL
#can now get values of specific species and sample eg. df["Aspergillus_flavus","ITS_iSeq"]



#calculate %ID and append result as new column to df
df[,"PE_ITS_iSeq_ID_rate"] <- (df[,"PE_ITS_iSeq_ID"]/(df[,"PE_ITS_iSeq_ID"] + df[,"PE_ITS_iSeq_Missing"]))*100

#calculate %ID and append result as new column to df
df[,"PE_ITS_iSeqrep2_ID_rate"] <- (df[,"PE_ITS_iSeqrep2_ID"]/(df[,"PE_ITS_iSeqrep2_ID"] + df[,"PE_ITS_iSeqrep2_Missing"]))*100

#calculate %ID and append result as new column to df
df[,"PE_ITS_MiSeq_ID_rate"] <- (df[,"PE_ITS_MiSeq_ID"]/(df[,"PE_ITS_MiSeq_ID"] + df[,"PE_ITS_MiSeq_Missing"]))*100

#TEF
df[,"merge_TEF_iSeq_IDrate"] <- (df[,"merge_TEF_iSeq_ID"]/(df[,"merge_TEF_iSeq_ID"] + df[,"merge_TEF_iSeq_Missing"]))*100

df[,"merge_TEF_MiSeq_IDrate"] <- (df[,"merge_TEF_MiSeq_ID"]/(df[,"merge_TEF_MiSeq_ID"] + df[,"merge_TEF_MiSeq_Missing"]))*100






write.csv(df,"ID_rate.csv", row.names=T) 