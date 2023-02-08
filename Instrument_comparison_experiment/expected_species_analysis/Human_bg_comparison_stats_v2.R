###### input expected only merged count files (ie. those with their proportions adjusted - performed for chi sq etc) with both TEF and ITS1 data
#only include mixes with low at 2%
file_list=list.files(pattern='.*merge.csv')
library(stringr)
library(tidyr)
library(dplyr)

####create matrix with expected proportion and PA values for all data ####

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

#### human vs non human ####################################################

#create column to indicate if human background present or not

human_rows <- grep("5.",all_data$Sample)

all_data[human_rows,6] <- "Yes"
all_data[-human_rows,6] <- "No"
colnames(all_data)[6] <- "Human_background"

###TEF######
#just TEF
TEF_all <- filter(all_data,grepl("TEF", Target))
range(TEF_all$Percent_abundance,na.rm=TRUE)

#convert species col from character to factor
TEF_all$Species <- as.factor(TEF_all$Species)
levels(TEF_all$Species)

#Assess data and reduce where necessary
#find those with at least 3 data points for each human bg status
TEF_Human_background_stats_summ <- TEF_all %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
Full_TEF_Human_background_stats_summ <- TEF_all %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance)
write.csv(Full_TEF_Human_background_stats_summ,file="full_TEF_stats_summ")


which(TEF_Human_background_stats_summ[,"n"]<3)
species_3ormore_ins <- TEF_Human_background_stats_summ[which(TEF_Human_background_stats_summ[,"n"]<3),"Species"]
TEF_r1 <- TEF_all[!TEF_all$Species %in% species_3ormore_ins$Species,]

#fudge to remove the levels of removed species
TEF_r1$Species <- as.character(TEF_r1$Species)
TEF_r1$Species <- as.factor(TEF_r1$Species)

rownames(TEF_r1) <- seq(length=nrow(TEF_r1))
TEF_r1$Human_background <- as.factor(TEF_r1$Human_background)
levels(TEF_r1$Human_background)

r1_summ <- TEF_r1 %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")




#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(r1_summ[,"n"] %% 2) == 0)

#Replace any 0 with NA
#rows_to_remove <- which(TEF_r1[,"Percent_abundance"] == 0)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  TEF_r1 <- TEF_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- TEF_r1 %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  TEF_r1 <- TEF_r1[!TEF_r1$Species %in% species_3ormore_ins$Species,]
}



TEF_r1$Species <- as.character(TEF_r1$Species)
TEF_r1$Species <- as.factor(TEF_r1$Species)

rownames(TEF_r1) <- seq(length=nrow(TEF_r1))
TEF_r1$Human_background <- as.factor(TEF_r1$Human_background)
levels(TEF_r1$Human_background)

cleaned_data_summary <- TEF_r1 %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")


#perform wilcox 
TEF_r1 %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Human_background,alt="two.sided")$p.value)
TEF_clean_wilcox <- TEF_r1 %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Human_background,alt="two.sided")$p.value)
write.csv(TEF_clean_wilcox,file="TEF_clean_wilcox_pvals")


#Plot cleaned data with wilcox results (using standard ggplot and stat_compare_means)
ggplot(TEF_r1,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Human_background)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Human\nbackground", guide = guide_legend(reverse=TRUE))+
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
ggplot(TEF_all,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Human_background)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Human\nbackground", guide = guide_legend(reverse=TRUE))+
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

####ITS1######
#just ITS1

ITS1_all <- filter(all_data,grepl("ITS1", Target))
range(ITS1_all$Percent_abundance,na.rm=TRUE)

#convert species col from character to factor
ITS1_all$Species <- as.factor(ITS1_all$Species)
levels(ITS1_all$Species)

###as I now have ITS iSeqrep2, need to split into subsections before applying cleaning & wilcox tests
#####split into iSeq & MiSeq####
ITS1_iSeqMiSeq <- filter(ITS1_all,!grepl("iSeqrep2",Instrument))


#Assess data and reduce where necessary
#find those with at least 3 data points for each human bg status
ITS1_iSeqMiSeq_Human_background_stats_summ <- ITS1_iSeqMiSeq %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

Full_ITS1_iSeqMiSeq_Human_background_stats_summ <- ITS1_iSeqMiSeq %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance)
write.csv(Full_ITS1_iSeqMiSeq_Human_background_stats_summ,file="full_ITS1_iSeqMiSeq_stats_summ")



which(ITS1_iSeqMiSeq_Human_background_stats_summ[,"n"]<3)
species_3ormore_ins <- ITS1_iSeqMiSeq_Human_background_stats_summ[which(ITS1_iSeqMiSeq_Human_background_stats_summ[,"n"]<3),"Species"]
ITS1_iSeqMiSeq_clean <- ITS1_iSeqMiSeq[!ITS1_iSeqMiSeq$Species %in% species_3ormore_ins$Species,]

#fudge to remove the levels of removed species
ITS1_iSeqMiSeq_clean$Species <- as.character(ITS1_iSeqMiSeq_clean$Species)
ITS1_iSeqMiSeq_clean$Species <- as.factor(ITS1_iSeqMiSeq_clean$Species)

rownames(ITS1_iSeqMiSeq_clean) <- seq(length=nrow(ITS1_iSeqMiSeq_clean))
ITS1_iSeqMiSeq_clean$Human_background <- as.factor(ITS1_iSeqMiSeq_clean$Human_background)
levels(ITS1_iSeqMiSeq_clean$Human_background)

ITS1_iSeqMiSeq_clean_summ <- ITS1_iSeqMiSeq_clean %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")


#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(ITS1_iSeqMiSeq_clean_summ[,"n"] %% 2) == 0)

#Replace any 0 with NA
#rows_to_remove <- which(ITS1_r1[,"Percent_abundance"] == 0)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  ITS1_r1 <- ITS1_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- ITS1_iSeqMiSeq_clean %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_iSeqMiSeq_clean <- ITS1_iSeqMiSeq_clean[!ITS1_iSeqMiSeq_clean$Species %in% species_3ormore_ins$Species,]
}



ITS1_iSeqMiSeq_clean$Species <- as.character(ITS1_iSeqMiSeq_clean$Species)
ITS1_iSeqMiSeq_clean$Species <- as.factor(ITS1_iSeqMiSeq_clean$Species)

rownames(ITS1_iSeqMiSeq_clean) <- seq(length=nrow(ITS1_iSeqMiSeq_clean))
ITS1_iSeqMiSeq_clean$Human_background <- as.factor(ITS1_iSeqMiSeq_clean$Human_background)
levels(ITS1_iSeqMiSeq_clean$Human_background)

iSeqMiSeq_cleaned_data_summary <- ITS1_iSeqMiSeq_clean %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")


#perform wilcox 
ITS1_iSeqMiSeq_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Human_background,alt="two.sided")$p.value)
ITS1_iMi_clean_wilcox <- ITS1_iSeqMiSeq_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Human_background,alt="two.sided")$p.value)
write.csv(ITS1_iMi_clean_wilcox,file="ITS1_iSeqMiSeq_clean_wilcox_pvals")


#Plot cleaned data with wilcox results (using standard ggplot and stat_compare_means)
tiff("Human_background_comparison_ITS1_iSeqMiSeq_clean_PA.tiff",width = 1000, height = 1000)
print(ggplot(ITS1_iSeqMiSeq_clean,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Human_background)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Human\nbackground", guide = guide_legend(reverse=TRUE))+
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


#####split into iSeq2 only ####
ITS1_iSeq2<- filter(ITS1_all,grepl("iSeqrep2",Instrument))


#Assess data and reduce where necessary
#find those with at least 3 data points for each human bg status
ITS1_iSeq2_Human_background_stats_summ <- ITS1_iSeq2 %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")

Full_ITS1_iSeq2_Human_background_stats_summ <- ITS1_iSeq2 %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance)
write.csv(Full_ITS1_iSeq2_Human_background_stats_summ,file="full_ITS1_iSeq2_stats_summ")



which(ITS1_iSeq2_Human_background_stats_summ[,"n"]<3)
species_3ormore_ins <- ITS1_iSeq2_Human_background_stats_summ[which(ITS1_iSeq2_Human_background_stats_summ[,"n"]<3),"Species"]
ITS1_iSeq2_clean <- ITS1_iSeq2[!ITS1_iSeq2$Species %in% species_3ormore_ins$Species,]

#fudge to remove the levels of removed species
ITS1_iSeq2_clean$Species <- as.character(ITS1_iSeq2_clean$Species)
ITS1_iSeq2_clean$Species <- as.factor(ITS1_iSeq2_clean$Species)

rownames(ITS1_iSeq2_clean) <- seq(length=nrow(ITS1_iSeq2_clean))
ITS1_iSeq2_clean$Human_background <- as.factor(ITS1_iSeq2_clean$Human_background)
levels(ITS1_iSeq2_clean$Human_background)

ITS1_iSeq2_clean_summ <- ITS1_iSeq2_clean %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")


#check for even number samples per species (if uneven, likely a mix missing somewhere..)
which(!(ITS1_iSeq2_clean_summ[,"n"] %% 2) == 0)

#Replace any 0 with NA
#rows_to_remove <- which(ITS1_r1[,"Percent_abundance"] == 0)
#if (!is_empty(rows_to_remove)) {
#remove all lines with NA - NB this means cannot do paired test, as some partner data will be missing
#  ITS1_r1 <- ITS1_r1[-rows_to_remove,]
#}


#still 3 or more reps each?
st_summ <- ITS1_iSeq2_clean %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")
#any have less than 3 reps after cleaning?
which(st_summ[,"n"]<3)
#if so, remove those species
if (!is_empty(which(st_summ[,"n"]<3))) {
  species_3ormore_ins <- st_summ[which(st_summ[,"n"]<3),"Species"]
  ITS1_iSeq2_clean <- ITS1_iSeq2_clean[!ITS1_iSeq2_clean$Species %in% species_3ormore_ins$Species,]
}



ITS1_iSeq2_clean$Species <- as.character(ITS1_iSeq2_clean$Species)
ITS1_iSeq2_clean$Species <- as.factor(ITS1_iSeq2_clean$Species)

rownames(ITS1_iSeq2_clean) <- seq(length=nrow(ITS1_iSeq2_clean))
ITS1_iSeq2_clean$Human_background <- as.factor(ITS1_iSeq2_clean$Human_background)
levels(ITS1_iSeq2_clean$Human_background)

iSeqMiSeq_cleaned_data_summary <- ITS1_iSeq2_clean %>% select(Species,Percent_abundance,Human_background) %>% group_by(Species, Human_background)  %>% get_summary_stats(Percent_abundance, type="mean_sd")


#perform wilcox 
ITS1_iSeq2_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Human_background,alt="two.sided")$p.value)
ITS1_iMi_clean_wilcox <- ITS1_iSeq2_clean %>%  group_by(Species) %>% summarise(p=wilcox.test(Percent_abundance~Human_background,alt="two.sided")$p.value)
write.csv(ITS1_iMi_clean_wilcox,file="ITS1_iSeq2_clean_wilcox_pvals")


#Plot cleaned data with wilcox results (using standard ggplot and stat_compare_means)
tiff("Human_background_comparison_ITS1_iSeq2_clean_PA.tiff",width = 1000, height = 1000)
print(ggplot(ITS1_iSeq2_clean,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Human_background)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Human\nbackground", guide = guide_legend(reverse=TRUE))+
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

####plot all ITS1 #### 
ggplot(ITS1_all,aes(x=gsub("_"," ",Species),y=Percent_abundance,fill=Human_background)) + 
  geom_boxplot(outlier.size=0.8,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.2, dotsize=0.3, position=position_dodge(1)) +
  coord_flip() +
  scale_fill_discrete(name="Human\nbackground", guide = guide_legend(reverse=TRUE))+
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


######### plot TEF alone, clean only ####################################


#perform wilcox just for TEF
TEF_r1$Species <- gsub("_"," ",TEF_r1$Species)
TEF.clean.stat <- compare_means(Percent_abundance ~ Human_background, data=TEF_r1, method= "wilcox.test",group.by=c("Species"))

#order species of full cleaned data
TEF_r1 <- TEF_r1[order(TEF_r1$Species,decreasing = TRUE),]


##plot using ggpubr and manually add stat info, with ND where appropriate
tiff("TEF_Human_background_clean_comparison_PA.tiff",width = 600, height = 800)
print(ggboxplot(TEF_r1, x = "Species", y = "Percent_abundance", width=0.7,
                fill = "Human_background", palette = c("#999999", "#E69F00"),
                add = "jitter", add.params = list(size = 0.4), sort.val = "desc") +
        #add stats manually
        stat_pvalue_manual(TEF.clean.stat, x="Species", y.position = 20,
                           label = "p.signif", position = position_dodge(0.8)) +
        #flip coords
        coord_flip() +
        #change legend order and title
        guides(fill = guide_legend(reverse = TRUE)) +
        labs(fill = "Human background") +
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

######## plot TEF alone - all ##########################################


#perform wilcox just for TEF
TEF_all$Species <- gsub("_"," ",TEF_all$Species)
###perform wilcox on all data, then replace significance of any with only 2 reps to 'ND'
TEF.stat <- compare_means(Percent_abundance ~ Human_background, data=TEF_all, method= "wilcox.test",group.by=c("Species"))



#find species with <3 reps 
TEF_sp_to_replace <- TEF_Human_background_stats_summ[which(TEF_Human_background_stats_summ[,"n"]<3),"Species"]
write.csv(TEF_sp_to_replace,file="TEF_sp_to_replace.csv")

##write to csv and manually edit all.stat
write.csv(TEF.stat,file="TEF_all_stats.csv")
#read in new modified stat table
TEF.stat.mod <- read.csv("TEF_all_stats_mod.csv",header=TRUE)

#order species of full, uncleaned data
TEF_all <- TEF_all[order(TEF_all$Species,decreasing = TRUE),]


##plot using ggpubr and manually add stat info, with ND where appropriate
tiff("TEF_Human_background_comparison_PA.tiff",width = 800, height = 1000)
print(ggboxplot(TEF_all, x = "Species", y = "Percent_abundance", width=0.7,
                fill = "Human_background", palette = c("#999999", "#E69F00"),
                add = "jitter", add.params = list(size = 0.4), sort.val = "desc") +
        #add stats manually
        stat_pvalue_manual(TEF.stat.mod, x="Species", y.position = 20,
                           label = "p.signif", position = position_dodge(0.8)) +
        #flip coords
        coord_flip() +
        #change legend order and title
        guides(fill = guide_legend(reverse = TRUE)) +
        labs(fill = "Human background") +
        #add lines
        geom_hline(yintercept = 2, colour = "grey40", size = 0.5,linetype="dashed") +
        geom_hline(yintercept = 0, colour = "red", size = 0.5,linetype="dotted") +
        #add lines between species
        geom_vline(xintercept=c(seq(1.5,16.5,by=1)),colour="grey",size=0.5)+
        xlab("Species") +
        ylab("Percent abundance")+
        theme(axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = .5, vjust = .5),
              axis.text.y = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = .5, face='italic'),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 25, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 25, angle = 90, hjust = .5, vjust = .5),
              panel.grid.minor.x=element_line(colour="grey", size=0.5),
              legend.title=element_text(size=22),
              legend.text=element_text(size=20),
              strip.background =element_rect(fill="grey"),
              strip.text=element_text(size=25, color='black')))
dev.off()

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
