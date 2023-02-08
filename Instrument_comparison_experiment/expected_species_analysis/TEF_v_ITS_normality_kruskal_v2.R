###### input expected only merged count files (ie. those with their proportions adjusted - performed for chi sq etc) with both TEF and ITS1 data
#only include mixes with low at 2%

file_list=list.files(pattern='.*merge.csv')
library(stringr)
library(tidyr)
library(dplyr)

#create matrix with expected proportion and PA values for all data

for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  MC=read.csv(filename, header=FALSE, sep=",",col.names=c("Counts","Species","Sample","Proportion"))
  MC <- MC[-1]
  by_sample <- spread(MC,Sample,Proportion)
  Mock_col <- grep(".*Mock", names(by_sample), value=T)
  
  mock_col_num <- which( colnames(by_sample)==Mock_col)
  colnames(by_sample)[mock_col_num] <- "Expected"
  #####Analyse PA for low species only
  ##create new data with only those species expected at low level (2%)
  rows_to_keep <- which(by_sample[,"Expected"]==0.02)
  two_percent_exp_only <- by_sample[rows_to_keep,]
  #remove expected proportions col
  two_percent_exp_only <- two_percent_exp_only[,-mock_col_num]
  ##make data long
  data_long <- two_percent_exp_only %>% gather(Sample,two_percent_exp_only,-Species)
  #convert propotions to percentage
  data_long[3] <- data_long[3]*100
  colnames(data_long)[3] <- "Percent_abundance"
  #write to file using samplename 
  samplename=str_replace(filename,"_merge.csv","")
  write.csv(data_long, paste(samplename,"_Percent_abundance_0.02_only.csv", sep=""), row.names=F)
}


#merge all data long mix PA files
filenames <- list.files(pattern='.*_Percent_abundance_0.02_only.csv')
tbl <- lapply(filenames, function(x){
  read.table(x, sep=",", as.is=TRUE, header=TRUE)})
all_data <- do.call(rbind, tbl)
#split sample column into mix name, instrument & target 
all_data <- separate(data = all_data, col = Sample, into = c("Sample", "Instrument","Target"), sep = "_") 
#replace ITS with ITS1
ITS_rows <- grep("ITS",all_data$Target)
#change ITS to ITS1 (for plotting labels)
all_data[ITS_rows,"Target"] <- "ITS1"


library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)

#############First, check normality & perform stats to check if PA varies between species (NOT INSTRUMENT COMPARISON)

##########################TEF############################
TEF_all <- filter(all_data,grepl("TEF", Target))

#convert species col from character to factor
TEF_all$Species <- as.factor(TEF_all$Species)
levels(TEF_all$Species)

TEF_all %>% select(Species,Percent_abundance) %>% group_by(Species) %>% get_summary_stats(Percent_abundance, type="mean_sd")
#Only take those with 3 or more data points
stats_summ <- TEF_all %>% select(Species,Percent_abundance) %>% group_by(Species) %>% get_summary_stats(Percent_abundance, type="mean_sd")
which(stats_summ[,"n"]<3)
species_3ormore <- stats_summ[which(stats_summ[,"n"]>=3),"Species"]

TEF_reduced <- TEF_all[TEF_all$Species %in% species_3ormore$Species,]
#fudge to remove the levels of removed species
TEF_reduced$Species <- as.character(TEF_reduced$Species)
TEF_reduced$Species <- as.factor(TEF_reduced$Species)
levels(TEF_reduced$Species)

ggboxplot(TEF_reduced,x="Species",y="Percent_abundance") +
  coord_flip()

TEF_reduced %>% select(Species,Percent_abundance) %>% group_by(Species) %>% identify_outliers()
#no extreme outliers

#######check normality assumptions 
##by analysing the model residuals
# Build the linear model
model  <- lm(Percent_abundance ~ Species, data = TEF_reduced)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
#p val was significant, so we cannot assume normality

########check normality by groups
shap_results <- TEF_reduced %>% select(Species,Percent_abundance) %>% group_by(Species) %>% shapiro_test(Percent_abundance)
#will error if sample size <3 for any species
#find sig p vals (ie. non normal species data)
shap_results[shap_results$p<0.05,]
#5 species not normal - the main 4 plus A flavus!

########as it doesnt seem normal, should use Kruskal-Wallis test
####Kruskal-Wallis test
#get more summary stats, including median, iqr etc
TEF_reduced %>% select(Species,Percent_abundance) %>% group_by(Species) %>% get_summary_stats(Percent_abundance, type="common")
#use the pipe-friendly kruskal_test() function [rstatix package], a wrapper around the R base function kruskal.test()
res.kruskal <- TEF_reduced %>% kruskal_test(Percent_abundance ~ Species)
res.kruskal
#perform Dunns for pairwise comparisons


###################################perform for ITS ##############################

#get just ITS1 data
ITS1_all <- filter(all_data,grepl("ITS1",Target))

#Only take those with 3 or more data points
stats_summ <- ITS1_all %>% select(Species,Percent_abundance) %>% group_by(Species) %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any with <3 reps?? (NB - all ITS, not for instrument comparison)
which(stats_summ[,"n"]<3)

species_3ormore <- stats_summ[which(stats_summ[,"n"]>=3),"Species"]

ITS1_reduced <- ITS1_all[ITS1_all$Species %in% species_3ormore$Species,]
#fudge to remove the levels of removed species
ITS1_reduced$Species <- as.character(ITS1_reduced$Species)
ITS1_reduced$Species <- as.factor(ITS1_reduced$Species)

ggboxplot(ITS1_reduced,x="Species",y="Percent_abundance") +
  coord_flip()

ITS1_reduced %>% select(Species,Percent_abundance) %>% group_by(Species) %>% identify_outliers()
#some extreme outliers

#####check normality assumptions 
##by analysing the model residuals
# Build the linear model
model  <- lm(Percent_abundance ~ Species, data = ITS1_reduced)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
#p val was significant, so we cannot assume normality

#will error if sample size <3 for any species
#######check normality by groups
shap_results <- ITS1_reduced %>% select(Species,Percent_abundance) %>% group_by(Species) %>% shapiro_test(Percent_abundance)
#will error if sample size <3 for any species
#find sig p vals (ie. non normal species data)
shap_results[shap_results$p<0.05,]


######as it doesnt seem normal, should use Kruskal-Wallis test
###Kruskal-Wallis test
#get more summary stats, including median, iqr etc
ITS1_reduced %>% select(Species,Percent_abundance) %>% group_by(Species) %>% get_summary_stats(Percent_abundance, type="common")
#use the pipe-friendly kruskal_test() function [rstatix package], a wrapper around the R base function kruskal.test()
res.ITS1.kruskal <- ITS1_reduced %>% kruskal_test(Percent_abundance ~ Species)
res.ITS1.kruskal$p

res.kruskal$p
