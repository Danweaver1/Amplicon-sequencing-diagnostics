## Performs chi squared tests on expected vs. experimental results for MiSeq data

file_list=list.files(pattern='.*Mock_proportional')
library(stringr)

for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  mock_counts <- read.table(filename,header=TRUE,sep=",")
  #get sample file names
  samplename=str_replace(filename,"_Mock_proportional.csv","")
  TEF_MiSeq=paste("TEF_Exp_sp_",samplename, "_MiSeq_cutoff_filtered_0_2_prop_adjusted.csv", sep="")
  ITS_MiSeq=paste("ITS_Exp_sp_",samplename, "_MiSeq_cutoff_filtered_0_2_prop_adjusted.csv", sep="")
  #check they both exist (rep1)
  if (file.exists(TEF_MiSeq)) {
    #read in files
    TEFsample_counts <- read.table(TEF_MiSeq,header=TRUE,sep=",")
    if (sum(TEFsample_counts[,4] != 0 )) {
      ####TEF chi with sample proportional counts
      TEF_res_pro <- chisq.test(TEFsample_counts[,4],p=mock_counts[,4])
      pvals<- c(TEF_res_pro$p.value)
      pval_table <- as.data.frame(pvals,row.names=c("TEFpro"))
      pval_table[2] <- rep(samplename,1)
      colnames(pval_table)<- c("pvals","sample")
      write.csv(pval_table,paste("TEF_",samplename,"_pvals.csv",sep=""), row.names=T) 
    } else {
      print(paste("Empty count files found for TEF",samplename,sep=" "))
    } 
  } else {
    print(paste("All sample files not available for TEF",samplename,sep=" "))
  }
  if (file.exists(ITS_MiSeq)) {
    #read in files
    ITSsample_counts <- read.table(ITS_MiSeq,header=TRUE,sep=",")
    if (sum(ITSsample_counts[,4] != 0)) {
      
      ###ITS chi with sample proportional counts
      ITS_res_pro <- chisq.test(ITSsample_counts[,4],p=mock_counts[,4])
      
      pvals<- c(ITS_res_pro$p.value)
      pval_table <- as.data.frame(pvals,row.names=c("ITSpro"))
      pval_table[2] <- rep(samplename,1)
      colnames(pval_table)<- c("pvals","sample")
      write.csv(pval_table,paste("ITS_", samplename,"_pvals.csv",sep=""), row.names=T) 
    } else {
      print(paste("Empty count files found for ITS",samplename,sep=" "))
    } 
  } else {
    print(paste("All sample files not available for ITS",samplename,sep=" "))
  }
}

temp = list.files(pattern="*_pvals.csv")
named.list <- lapply(temp, read.csv)
library(data.table)
files.matrix <-rbindlist(named.list)
colnames(files.matrix) = c("test","pvals","sample")
write.csv(files.matrix,"matrix_p_vals.csv", row.names=F) 

