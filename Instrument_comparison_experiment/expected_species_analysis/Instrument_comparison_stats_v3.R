###### input expected only merged count files (ie. those with their proportions adjusted - performed for chi sq etc) with both TEF and ITS1 data
#### setup####
#only include mixes with low at 2%
file_list=list.files(pattern='.*merge.csv')
library(stringr)
library(tidyr)
library(dplyr)
 
#######create matrix with expected proportion and PA values for all data####

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


#######merge all data long mix PA files####
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



#######compare instrument for each species ##########################################

######First, need to check if, for each target, any species have <3 reps per instrument#######

#####TEF####
#just TEF
TEF_all <- filter(all_data,grepl("TEF", Target))
range(TEF_all$Percent_abundance,na.rm=TRUE)
median(TEF_all$Percent_abundance)
#convert species col from character to factor
TEF_all$Species <- as.factor(TEF_all$Species)
levels(TEF_all$Species)

#Assess data and reduce where necessary
#find those with at least 3 data points for each instrument 
TEF_instrument_stats_summ <- TEF_all %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

#get instrument level summary stats
full_TEF_instrument_stats_summ <- TEF_all %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance)
write.csv(full_TEF_instrument_stats_summ,file="Full_TEF_instrument_stats_summ")
#get summary stats for both instruments combined
full_TEF_stats_summ <- TEF_all %>% select(Species,Percent_abundance) %>% group_by(Species)  %>% get_summary_stats(Percent_abundance)
write.csv(full_TEF_stats_summ,file="Full_TEF_stats_summ")


which(TEF_instrument_stats_summ[,"n"]<3)
species_3ormore_ins <- TEF_instrument_stats_summ[which(TEF_instrument_stats_summ[,"n"]<3),"Species"]

TEF_r1 <- TEF_all[!TEF_all$Species %in% species_3ormore_ins$Species,]
r1_summ <- TEF_r1 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

#fudge to remove the levels of removed species
TEF_r1$Species <- as.character(TEF_r1$Species)
TEF_r1$Species <- as.factor(TEF_r1$Species)

rownames(TEF_r1) <- seq(length=nrow(TEF_r1))

TEF_r1$Instrument <- as.factor(TEF_r1$Instrument)
levels(TEF_r1$Instrument)
#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(r1_summ[,"n"] %% 2) == 0)

#Replace any 0 with NA
#rows_to_remove <- which(TEF_r1[,"Percent_abundance"] == 0)
#if (!is_empty(rows_to_remove)) {
  #remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  TEF_r1 <- TEF_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- TEF_r1 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_r1 <- ITS1_r1[!ITS1_r1$Species %in% species_3ormore_ins$Species,]
}

TEF_cleaned_data_summary <- TEF_r1 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")


rownames(TEF_r1) <- seq(length=nrow(TEF_r1))

TEF_r1$Instrument <- as.factor(TEF_r1$Instrument)
levels(TEF_r1$Instrument)

#perform wilcox 
TEF_r1 %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
TEF_clean_wilcox <- TEF_r1 %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
write.csv(TEF_clean_wilcox,file="TEF_clean_wilcox_pvals")
#Plot cleaned data with wilcox results

ggplot(TEF_r1,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
        geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
        geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
        coord_flip() +
        scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
        geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
        geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
        #add stat sig
        stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
        xlab("Species") +
        ylab("Percent abundance") +
        #theme_light() +
        theme_linedraw() +
        theme(axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
              axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25),
              strip.background =element_rect(fill="grey"),
              strip.text=element_text(size=30, color='black'))

###plot all TEF 
ggplot(TEF_all,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
  geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
  geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
  #add stat sig
  stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
  xlab("Species") +
  ylab("Percent abundance") +
  #theme_light() +
  theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),
        strip.background =element_rect(fill="grey"),
        strip.text=element_text(size=30, color='black'))



#####ITS1####
####get just ITS1 data####
ITS1_all <- filter(all_data,grepl("ITS1",Target))
range(ITS1_all$Percent_abundance,na.rm=TRUE)
median(ITS1_all$Percent_abundance)

#convert species col from character to factor
ITS1_all$Species <- as.factor(ITS1_all$Species)
levels(ITS1_all$Species)

###as I now have ITS iSeqrep2, need to split into subsections before applying wilcox tests
####split into iSeq vs MiSeq####
ITS1_iSeqMiSeq <- filter(ITS1_all,!grepl("iSeqrep2",Instrument))

#Assess data and reduce where necessary
#find those with at least 3 data points for each instrument 
ITS_iSeqMiSeq_instrument_stats_summ <- ITS1_iSeqMiSeq %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#get instrument level summary stats
full_ITS_iSeqMiSeq_instrument_stats_summ <- ITS1_iSeqMiSeq %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_iSeqMiSeq_instrument_stats_summ,file="full_ITS_iSeqMiSeq_instrument_summ_stats")
#get summary stats for both instruments combined
full_ITS_iSeqMiSeq_stats_summ <- ITS1_iSeqMiSeq %>% select(Species,Percent_abundance) %>% group_by(Species)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_iSeqMiSeq_stats_summ,file="full_ITS_iSeqMiSeq_summ_stats")


which(ITS_iSeqMiSeq_instrument_stats_summ[,"n"]<3)
species_3ormore_ins <- ITS_iSeqMiSeq_instrument_stats_summ[which(ITS_iSeqMiSeq_instrument_stats_summ[,"n"]<3),"Species"]

ITS1_iSeqMiSeq_clean <- ITS1_iSeqMiSeq[!ITS1_iSeqMiSeq$Species %in% species_3ormore_ins$Species,]
ITS1_iSeqMiSeq_clean_summ <- ITS1_iSeqMiSeq_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

#fudge to remove the levels of removed species
ITS1_iSeqMiSeq_clean$Species <- as.character(ITS1_iSeqMiSeq_clean$Species)
ITS1_iSeqMiSeq_clean$Species <- as.factor(ITS1_iSeqMiSeq_clean$Species)
levels(ITS1_iSeqMiSeq_clean$Species)
#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(ITS1_iSeqMiSeq_clean_summ[,"n"] %% 2) == 0)

#Remove any 0 
#rows_to_remove <- which(ITS1_r1[,"Percent_abundance"] == 0)
#length(rows_to_remove)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  ITS1_r1 <- ITS1_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- ITS1_iSeqMiSeq_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_iSeqMiSeq_clean <- ITS1_iSeqMiSeq_clean[!ITS1_iSeqMiSeq_clean$Species %in% species_3ormore_ins$Species,]
}

ITS_iSeqMiSeq_cleaned_data_summary <- ITS1_iSeqMiSeq_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

rownames(ITS1_iSeqMiSeq_clean) <- seq(length=nrow(ITS1_iSeqMiSeq_clean))

ITS1_iSeqMiSeq_clean$Instrument <- as.factor(ITS1_iSeqMiSeq_clean$Instrument)
levels(ITS1_iSeqMiSeq_clean$Instrument)





#perform wilcox 
ITS1_iSeqMiSeq_clean %>% group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
ITS1_iSeqMiSeq_clean_wilcox <- ITS1_iSeqMiSeq_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
write.csv(ITS1_iSeqMiSeq_clean_wilcox,file="ITS1_iSeqMiSeq_clean_wilcox_pvals")
#Plot cleaned data with wilcox results

x11()
tiff("instrument_comparison_clean_ITS_iSeqMiSeq_PA.tiff",width = 1000, height = 800)
print(ggplot(ITS1_iSeqMiSeq_clean,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
  geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
  geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
  #add stat sig
  stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
  xlab("Species") +
  ylab("Percent abundance") +
  #theme_light() +
  theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),
        strip.background =element_rect(fill="grey"),
        strip.text=element_text(size=30, color='black')))
dev.off()


####split into iSeq vs iSeq2####
ITS1_iSeqiSeq2 <- filter(ITS1_all,!grepl("MiSeq",Instrument))

#Assess data and reduce where necessary
#find those with at least 3 data points for each instrument 
ITS_iSeqiSeq2_instrument_stats_summ <- ITS1_iSeqiSeq2 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#get instrument level summary stats
full_ITS_iSeqiSeq2_instrument_stats_summ <- ITS1_iSeqiSeq2 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_iSeqiSeq2_instrument_stats_summ,file="full_ITS_iSeqiSeq2_instrument_summ_stats")
#get summary stats for both instruments combined
full_ITS_iSeqiSeq2_stats_summ <- ITS1_iSeqiSeq2 %>% select(Species,Percent_abundance) %>% group_by(Species)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_iSeqiSeq2_stats_summ,file="full_ITS_iSeqiSeq2_summ_stats")


which(ITS_iSeqiSeq2_instrument_stats_summ[,"n"]<3)
species_3ormore_ins <- ITS_iSeqiSeq2_instrument_stats_summ[which(ITS_iSeqiSeq2_instrument_stats_summ[,"n"]<3),"Species"]

ITS1_iSeqiSeq2_clean <- ITS1_iSeqiSeq2[!ITS1_iSeqiSeq2$Species %in% species_3ormore_ins$Species,]
ITS1_iSeqiSeq2_clean_summ <- ITS1_iSeqiSeq2_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

#fudge to remove the levels of removed species
ITS1_iSeqiSeq2_clean$Species <- as.character(ITS1_iSeqiSeq2_clean$Species)
ITS1_iSeqiSeq2_clean$Species <- as.factor(ITS1_iSeqiSeq2_clean$Species)
levels(ITS1_iSeqiSeq2_clean$Species)
#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(ITS1_iSeqiSeq2_clean_summ[,"n"] %% 2) == 0)

#Remove any 0 
#rows_to_remove <- which(ITS1_r1[,"Percent_abundance"] == 0)
#length(rows_to_remove)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  ITS1_r1 <- ITS1_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- ITS1_iSeqiSeq2_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_iSeqiSeq2_clean <- ITS1_iSeqiSeq2_clean[!ITS1_iSeqiSeq2_clean$Species %in% species_3ormore_ins$Species,]
}

ITS_iSeqiSeq2_cleaned_data_summary <- ITS1_iSeqiSeq2_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

rownames(ITS1_iSeqiSeq2_clean) <- seq(length=nrow(ITS1_iSeqiSeq2_clean))

ITS1_iSeqiSeq2_clean$Instrument <- as.factor(ITS1_iSeqiSeq2_clean$Instrument)
levels(ITS1_iSeqiSeq2_clean$Instrument)


#perform wilcox 
ITS1_iSeqiSeq2_clean %>% group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
ITS1_iSeqiSeq2_clean_wilcox <- ITS1_iSeqiSeq2_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
write.csv(ITS1_iSeqiSeq2_clean_wilcox,file="ITS1_iSeqiSeq2_clean_wilcox_pvals")
#Plot cleaned data with wilcox results

x11()
tiff("instrument_comparison_clean_ITS_iSeqiSeq2_PA.tiff",width = 1000, height = 800)
print(ggplot(ITS1_iSeqiSeq2_clean,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
  geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
  geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
  #add stat sig
  stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
  xlab("Species") +
  ylab("Percent abundance") +
  #theme_light() +
  theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),
        strip.background =element_rect(fill="grey"),
        strip.text=element_text(size=30, color='black')))
dev.off()

####split into iSeq2 vs MiSeq####
ITS1_iSeq2MiSeq <- filter(ITS1_all,grepl("iSeqrep2",Instrument) | grepl("MiSeq",Instrument) )

#Assess data and reduce where necessary
#find those with at least 3 data points for each instrument 
ITS_iSeq2MiSeq_instrument_stats_summ <- ITS1_iSeq2MiSeq %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#get instrument level summary stats
full_ITS_iSeq2MiSeq_instrument_stats_summ <- ITS1_iSeq2MiSeq %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_iSeq2MiSeq_instrument_stats_summ,file="full_ITS_iSeq2MiSeq_instrument_summ_stats")
#get summary stats for both instruments combined
full_ITS_iSeq2MiSeq_stats_summ <- ITS1_iSeq2MiSeq %>% select(Species,Percent_abundance) %>% group_by(Species)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_iSeq2MiSeq_stats_summ,file="full_ITS_iSeq2MiSeq_summ_stats")


which(ITS_iSeq2MiSeq_instrument_stats_summ[,"n"]<3)
species_3ormore_ins <- ITS_iSeq2MiSeq_instrument_stats_summ[which(ITS_iSeq2MiSeq_instrument_stats_summ[,"n"]<3),"Species"]

ITS1_iSeq2MiSeq_clean <- ITS1_iSeq2MiSeq[!ITS1_iSeq2MiSeq$Species %in% species_3ormore_ins$Species,]
ITS1_iSeq2MiSeq_clean_summ <- ITS1_iSeq2MiSeq_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

#fudge to remove the levels of removed species
ITS1_iSeq2MiSeq_clean$Species <- as.character(ITS1_iSeq2MiSeq_clean$Species)
ITS1_iSeq2MiSeq_clean$Species <- as.factor(ITS1_iSeq2MiSeq_clean$Species)
levels(ITS1_iSeq2MiSeq_clean$Species)
#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(ITS1_iSeq2MiSeq_clean_summ[,"n"] %% 2) == 0)

#Remove any 0 
#rows_to_remove <- which(ITS1_r1[,"Percent_abundance"] == 0)
#length(rows_to_remove)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  ITS1_r1 <- ITS1_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- ITS1_iSeq2MiSeq_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_iSeq2MiSeq_clean <- ITS1_iSeq2MiSeq_clean[!ITS1_iSeq2MiSeq_clean$Species %in% species_3ormore_ins$Species,]
}

ITS_iSeq2MiSeq_cleaned_data_summary <- ITS1_iSeq2MiSeq_clean %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

rownames(ITS1_iSeq2MiSeq_clean) <- seq(length=nrow(ITS1_iSeq2MiSeq_clean))

ITS1_iSeq2MiSeq_clean$Instrument <- as.factor(ITS1_iSeq2MiSeq_clean$Instrument)
levels(ITS1_iSeq2MiSeq_clean$Instrument)

#perform wilcox 
ITS1_iSeq2MiSeq_clean %>% group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
ITS1_iSeq2MiSeq_clean_wilcox <- ITS1_iSeq2MiSeq_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Instrument,alt="two.sided")$p.value)
write.csv(ITS1_iSeq2MiSeq_clean_wilcox,file="ITS1_clean_wilcox_pvals")
#Plot cleaned data with wilcox results
tiff("instrument_comparison_clean_ITS_iSeq2MiSeq_PA.tiff",width = 1000, height = 800)
print(ggplot(ITS1_iSeq2MiSeq_clean,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
  geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
  geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
  #add stat sig
  stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
  xlab("Species") +
  ylab("Percent abundance") +
  #theme_light() +
  theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),
        strip.background =element_rect(fill="grey"),
        strip.text=element_text(size=30, color='black')))
dev.off()

####plot all ITS#### 
x11()
tiff("instrument_comparison_all_ITS_PA.tiff",width = 1000, height = 800)
print(ggplot(ITS1_all,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
  geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
  geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
  #add stat sig
  stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
  xlab("Species") +
  ylab("Percent abundance") +
  #theme_light() +
  theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),
        strip.background =element_rect(fill="grey"),
        strip.text=element_text(size=30, color='black')))
dev.off()

####clean all ITS together ####

#Assess data and reduce where necessary
#find those with at least 3 data points for each instrument 
ITS_instrument_stats_summ <- ITS1_all %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#get instrument level summary stats
full_ITS_instrument_stats_summ <- ITS1_all %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_instrument_stats_summ,file="full_ITS_instrument_summ_stats")
#get summary stats for both instruments combined
full_ITS_stats_summ <- ITS1_all %>% select(Species,Percent_abundance) %>% group_by(Species)  %>% get_summary_stats(Percent_abundance)
write.csv(full_ITS_stats_summ,file="full_ITS_summ_stats")


which(ITS_instrument_stats_summ[,"n"]<3)
species_3ormore_ins <- ITS_instrument_stats_summ[which(ITS_instrument_stats_summ[,"n"]<3),"Species"]

ITS1_r1 <- ITS1_all[!ITS1_all$Species %in% species_3ormore_ins$Species,]
ITS1_r1_summ <- ITS1_r1 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

#fudge to remove the levels of removed species
ITS1_r1$Species <- as.character(ITS1_r1$Species)
ITS1_r1$Species <- as.factor(ITS1_r1$Species)
levels(ITS1_r1$Species)
#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(ITS1_r1_summ[,"n"] %% 2) == 0)

#Remove any 0 
#rows_to_remove <- which(ITS1_r1[,"Percent_abundance"] == 0)
#length(rows_to_remove)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  ITS1_r1 <- ITS1_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- ITS1_r1 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_r1 <- ITS1_r1[!ITS1_r1$Species %in% species_3ormore_ins$Species,]
}

ITS_iSeq2MiSeq_cleaned_data_summary <- ITS1_r1 %>% select(Species,Percent_abundance,Instrument) %>% group_by(Species, Instrument)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

rownames(ITS1_r1) <- seq(length=nrow(ITS1_r1))

ITS1_r1$Instrument <- as.factor(ITS1_r1$Instrument)
levels(ITS1_r1$Instrument)


#####merge all cleaned data & plot ##########################

full_clean <- rbind(TEF_r1,ITS1_r1)
#plot with ggplot standard method
ggplot(full_clean,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Instrument)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Sequencing\ninstrument", guide = guide_legend(reverse=TRUE))+
  geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
  geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
  facet_wrap(~Target) +
  #add stat sig
  stat_compare_means(method= "wilcox.test",paired=FALSE,aes(label = ..p.signif..),label.x = 1.5, label.y = 82)+
  xlab("Species") +
  ylab("Percent abundance") +
  #theme_light() +
  theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),
        strip.background =element_rect(fill="grey"),
        strip.text=element_text(size=30, color='black'))


#####plotting full data but with stats for only 3 rep or more species, giving others 'ND' ##########

#TEF & ITS1 cleaned data wilcox together - this is to check the stats which I want shown at end!
full_clean$Species <- gsub("_"," ",full_clean$Species)
full.stat <- compare_means(Percent_abundance ~ Instrument, data=full_clean, method= "wilcox.test",group.by=c("Species","Target"))

#replace underscore with space
all_data$Species <- gsub("_"," ",all_data$Species)
###perform wilcox on all data, then replace significance of any with only 2 reps to 'ND'
all.stat <- compare_means(Percent_abundance ~ Instrument, data=all_data, method= "wilcox.test",group.by=c("Species","Target"))

#find species with <3 reps 
TEF_sp_to_replace <- TEF_instrument_stats_summ[which(TEF_instrument_stats_summ[,"n"]<3),"Species"]
write.csv(TEF_sp_to_replace,file="TEF_sp_to_replace.csv")
ITS_sp_to_replace <- ITS_instrument_stats_summ[which(ITS_instrument_stats_summ[,"n"]<3),"Species"]
write.csv(ITS_sp_to_replace,file="ITS_sp_to_replace.csv")

##write to csv and manually edit all.stat
write.csv(all.stat,file="all_data_stats.csv")
#read in new modified stat table
mod.all.stat <- read.csv("all_data_stats_mod1.csv",header=TRUE)
#order species of full, uncleaned data
all_data <- all_data[order(all_data$Species,decreasing = TRUE),]

##plot using ggpubr and manually add stat info, with ND where appropriate
tiff("instrument_comparison_PA.tiff",width = 1000, height = 1000)
print(ggboxplot(all_data, x = "Species", y = "Percent_abundance", width=0.6,
          fill = "Instrument", palette = c("#ccebc5", "#1b9e77","#4597B0"),
          add = "jitter", facet.by = "Target", add.params = list(size = 0.4), sort.val = "desc") +
        #add stats manually
        stat_pvalue_manual(mod.all.stat, x="Species", y.position = 84,
         label = "p.signif", position = position_dodge(0.8)) +
        #flip coords
        coord_flip() +
        #change legend order and title
        guides(fill = guide_legend(reverse = TRUE)) +
        labs(fill = "Sequencing instrument") +
        #add lines
        geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
        geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
        #add lines between species
        geom_vline(xintercept=c(seq(1.5,16.5,by=1)),colour="grey",size=0.5)+
        xlab("Species") +
        ylab("Percent abundance")+
        theme(axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 22, angle =0, hjust = .5, vjust = .5),
              axis.text.y = element_text(color = "grey20", size = 22, angle =0, hjust = 1, vjust = .5, face='italic'),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 25, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 25, angle = 90, hjust = .5, vjust = .5),
              panel.grid.minor.x=element_line(colour="grey", size=0.5),
              legend.title=element_text(size=25),
              legend.text=element_text(size=22),
              strip.background =element_rect(fill="grey"),
              strip.text=element_text(size=25, color='black')))
dev.off()


#####plot only clean data as above ###############################
#TEF & ITS1 cleaned data wilcox together - this is to check the stats which I want shown at end!
full_clean$Species <- gsub("_"," ",full_clean$Species)
full.stat <- compare_means(Percent_abundance ~ Instrument, data=full_clean, method= "wilcox.test",group.by=c("Species","Target"))

#order species of full cleaned data
full_clean <- full_clean[order(full_clean$Species,decreasing = TRUE),]

##plot using ggpubr and manually add stat info, with ND where appropriate
tiff("instrument_comparison_clean_PA.tiff",width = 1000, height = 800)
print(ggboxplot(full_clean, x = "Species", y = "Percent_abundance", width=0.6,
                fill = "Instrument", palette = c("#ccebc5", "#1b9e77","#4597B0"),
                add = "jitter", facet.by = "Target", add.params = list(size = 0.4), sort.val = "desc") +
        #add stats manually
        stat_pvalue_manual(full.stat, x="Species", y.position = 84,
                           label = "p.signif", position = position_dodge(0.8)) +
        #flip coords
        coord_flip() +
        #change legend order and title
        guides(fill = guide_legend(reverse = TRUE)) +
        labs(fill = "Sequencing instrument") +
        #add lines
        geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
        geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
        #add lines between species
        geom_vline(xintercept=c(seq(1.5,16.5,by=1)),colour="grey",size=0.5)+
        xlab("Species") +
        ylab("Percent abundance")+
        theme(axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
              axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face='italic'),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              panel.grid.minor.x=element_line(colour="grey", size=0.5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25),
              strip.background =element_rect(fill="grey"),
              strip.text=element_text(size=30, color='black')))
dev.off()



#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
