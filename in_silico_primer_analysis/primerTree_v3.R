## Performs in silico analysis of primer hits using primerTree package ####
##  A restriction is placed on organisms searched - set by "organism" parameter - set to Fungi here
## NOTE! Requires an NCBI API key to run

#### setup ####
#install.packages("primerTree")
library(primerTree)
library(ape)
library(ggplot2)
library(dplyr)
library(rstatix)

#### ITS1 degenerate primer set, 2000 aligns, 50 permutations ####
###Run primer search
ITS1F <- search_primer_pair(name='ITS1Ftest','TCCGTAGGTGAACCTGCGG', 'GGCYRCGTTCTTCATCGAYGC',
                            total_primer_specificity_mismatch=0,num_aligns=10,num_permutations=8,
                            organism="Fungi", api_key="")


ITS1F$BLAST_result[2]
accession2gi(ITS1F$BLAST_result[3])


tiff(filename="ITS1F",height=1000,width=1000)
print(plot(ITS1F))
dev.off()

ITS1F_lengths <- seq_lengths(ITS1F,summarize=TRUE)

tiff(filename="ITS1F_lengths",height=1000,width=1000)
print(barplot(ITS1F_lengths,main = "Frequency of sequence lengths for ITS1 primers",xlab="Amplicon length (in bp)",ylab=("Frequency")))
dev.off()

ITS1F_Blast <- ITS1F$BLAST_result
write.csv(ITS1F_Blast,file="ITS1F_blast_results")

length(ITS1F_Blast$gi)
length(unique(ITS1F_Blast$gi))

median(ITS1F_Blast$product_length)
mean(ITS1F_Blast$product_length)
range(ITS1F_Blast$product_length)

ITS1F_summ <- get_summary_stats(ITS1F_Blast)
write.csv(ITS1F_summ,file="ITS1F_results_summary")


write.csv(ITS1F$taxonomy,file="ITS1F_taxonomy")
#unique classes
ITS_taxa <- ITS1F$taxonomy
length(unique(ITS_taxa$class))
uniq_ITS_taxa <- sort(unique(ITS_taxa$class))
write.csv(uniq_ITS_taxa,file="unique_ITS_classes")

write.tree(ITS1F$tree,file="ITS1F_tree")

##filter seqs for plotting (any above 0.8kb)
filtered_lengths_ITS <- ITS1F_Blast$product_length[ITS1F_Blast$product_length < 800]

range(filtered_lengths_ITS)
mean(filtered_lengths_ITS)
median(filtered_lengths_ITS)
write.csv(filtered_lengths_ITS,file="filtered_lengths_ITS")

tiff(filename="ITS1F_lengths_filtered",height=1000,width=1000)
hist(filtered_lengths_ITS,breaks=seq(from=0,to=500,by=10),main = "Frequency of sequence lengths for ITS1 primers",xlab="Amplicon length (in bp)",ylab="Frequency",xlim=c(0,500))
dev.off()

##plot trees ## 

#plot&print with legend on bottom
ITSp_leg <- plot_tree(ITS1F$tree,type="fan",rank="class",taxonomy=ITS1F$taxonomy, size = 4,legend_cutoff=50,main = "ITS1")
ITSp_leg <- ITSp_leg +theme_void() + theme(text = element_text(size=30),legend.position = "bottom")

tiff(filename="ITS1F_tree_legend1",height=1000,width=1600)
print(ITSp_leg)
dev.off()

  #plot&print with legend on right
ITSp_leg2 <- plot_tree(ITS1F$tree,type="fan",rank="class",taxonomy=ITS1F$taxonomy, size = 4,legend_cutoff =50,main = "ITS1")
ITSp_leg2 <- ITSp_leg2 +theme_void() + theme(text = element_text(size=30))


tiff(filename="ITS1F_tree_legend2",height=600,width=800)
print(ITSp_leg2)
dev.off()

#plot&print without legend
ITSp <- plot_tree(ITS1F$tree,type="fan",rank="class",taxonomy=ITS1F$taxonomy, size = 1,legend_cutoff =50,main = "ITS1")
ITSp <- ITSp +theme_void() + theme(legend.position = "none")


tiff(filename="ITS1F_tree.tiff",height=650,width=600)
print(ITSp)
dev.off()

#plot&print with legend on bottom
ITSp_leg <- plot_tree(ITS1F$tree,type="fan",rank="phylum",taxonomy=ITS1F$taxonomy, size = 4,legend_cutoff=50,main = "ITS1")
ITSp_leg <- ITSp_leg +theme_void() + theme(text = element_text(size=30),legend.position = "bottom")

tiff(filename="ITS1F_tree_phylum_legend1",height=1000,width=1600)
print(ITSp_leg)
dev.off()

#plot&print with legend on right
ITSp_leg2 <- plot_tree(ITS1F$tree,type="fan",rank="phylum",taxonomy=ITS1F$taxonomy, size = 4,legend_cutoff =50,main = "ITS1")
ITSp_leg2 <- ITSp_leg2 +theme_void() + theme(text = element_text(size=30))


tiff(filename="ITS1F_tree_phylum_legend2",height=600,width=800)
print(ITSp_leg2)
dev.off()

#plot&print without legend
ITSp <- plot_tree(ITS1F$tree,type="fan",rank="phylum",taxonomy=ITS1F$taxonomy, size = 1,legend_cutoff =50,main = "ITS1")
ITSp <- ITSp +theme_void() + theme(legend.position = "none")


tiff(filename="ITS1F_tree_phylum.tiff",height=650,width=600)
print(ITSp)
dev.off()


#plot&print with legend on bottom
ITSp_leg <- plot_tree(ITS1F$tree,type="fan",rank="phylum",taxonomy=ITS1F$taxonomy, size = 4,legend_cutoff=50,main = "ITS1")
ITSp_leg <- ITSp_leg +theme_void() + theme(text = element_text(size=30),legend.position = "bottom")

tiff(filename="ITS1F_tree_fan_phylum_legend1",height=1000,width=1600)
print(ITSp_leg)
dev.off()

#plot&print with legend on right
ITSp_leg2 <- plot_tree(ITS1F$tree,type="fan",rank="phylum",taxonomy=ITS1F$taxonomy, size = 4,legend_cutoff =50,main = "ITS1")
ITSp_leg2 <- ITSp_leg2 +theme_void() + theme(text = element_text(size=30))


tiff(filename="ITS1F_tree_fan_phylum_legend2",height=600,width=800)
print(ITSp_leg2)
dev.off()

#plot&print without legend
ITSp <- plot_tree(ITS1F$tree,type="fan",rank="phylum",taxonomy=ITS1F$taxonomy, size = 1,legend_cutoff =50,main = "ITS1")
ITSp <- ITSp +theme_void() + theme(legend.position = "none")


tiff(filename="ITS1F_tree_phylum.tiff",height=650,width=600)
print(ITSp)
dev.off()
#plot&print fan tree without legend or labels
ITSp2 <- plot_tree(ITS1F$tree,type="unrooted",rank="phylum",taxonomy=ITS1F$taxonomy, size = 2.5,legend_cutoff =0,main = "ITS1") 
ITSp2 <- ITSp2 +theme_void() + theme(legend.position = "none")


tiff(filename="ITS1F_unrooted_tree_phylum",height=800,width=800)
print(ITSp2)
dev.off()

  #plot&print phylogram without legend
ITSp3 <- plot_tree(ITS1F$tree,type="phylogram",rank="class",taxonomy=ITS1F$taxonomy, size = 1,legend_cutoff =50,main = "ITS1")
ITSp3 <- ITSp3 +theme_void() + theme(legend.position = "none")


tiff(filename="ITS1F_phylo_tree",height=200,width=350)
print(ITSp3)
dev.off()

####TEF merge  ####

TEF_mergeF <- search_primer_pair(name='TEF_mergeF','YGGYGARTTCGARGCYGGTAT', 'ADCATRTTRTCDCCRTKGMADCC',
                                 total_primer_specificity_mismatch=3,num_aligns=1000,
                                 num_permutations=50,organism="Fungi", 
                                 api_key="")


tiff(filename="TEF_mergeF",height=1000,width=1000)
print(plot(TEF_mergeF))
dev.off()

TEF_mergeF_lengths <- seq_lengths(TEF_mergeF,summarize=TRUE)

TEF_mergeF_Blast <- TEF_mergeF$BLAST_result
write.csv(TEF_mergeF_Blast,file="TEF_mergeF_blast_results")

length(TEF_mergeF_Blast$gi)
length(unique(TEF_mergeF_Blast$gi))

median(TEF_mergeF_Blast$product_length)
mean(TEF_mergeF_Blast$product_length)
range(TEF_mergeF_Blast$product_length)

TEF_mergeF_summ <- get_summary_stats(TEF_mergeF_Blast)
write.csv(TEF_mergeF_summ,file="TEF_mergeF_results_summary")


write.csv(TEF_mergeF$taxonomy,file="TEF_mergeF_taxonomy")

#unique classes
TEF_taxa <- TEF_mergeF$taxonomy
length(unique(TEF_taxa$class))
uniq_TEF_taxa <- sort(unique(TEF_taxa$class))
write.csv(uniq_TEF_taxa,file="unique_TEF_classes")

write.tree(TEF_mergeF$tree,file="TEF_mergeF_tree")

tiff(filename="TEF_mergeF_lengths2",height=600,width=800)
hist(TEF_mergeF_Blast$product_length,breaks=seq(from=50,to=850,by=10),main = "Frequency of sequence lengths for TEF primers",xlab="Amplicon length (in bp)",ylab="Frequency",xlim=c(0,900))
dev.off()

##filter seqs for plotting (any above 0.8kb)
filtered_lengths_TEF <- TEF_mergeF_Blast$product_length[TEF_mergeF_Blast$product_length < 800]

range(filtered_lengths_TEF)
median(filtered_lengths_TEF)
mean(filtered_lengths_TEF)
write.csv(filtered_lengths_TEF,file="filtered_lengths_TEF")

tiff(filename="TEF_lengths_filtered",height=1000,width=1000)
hist(filtered_lengths_TEF,breaks=seq(from=0,to=500,by=10),main = "Frequency of sequence lengths for TEF primers",xlab="Amplicon length (in bp)",ylab="Frequency",xlim=c(0,500))
dev.off()

##plot trees ##

#plot&print with legend on bottom
TEFp_leg <- plot_tree(TEF_mergeF$tree,type="fan",rank="class",taxonomy=TEF_mergeF$taxonomy, size = 4,legend_cutoff =50,main = "TEF")
TEFp_leg <- TEFp_leg +theme_void() + theme(text = element_text(size=30),legend.position = "bottom")

tiff(filename="TEF_mergeF_tree_legend1",height=1000,width=1600)
print(TEFp_leg)
dev.off()

#plot&print with legend on right
TEFp_leg2 <- plot_tree(TEF_mergeF$tree,type="fan",rank="class",taxonomy=TEF_mergeF$taxonomy, size = 4,legend_cutoff =50,main = "TEF")
TEFp_leg2 <- TEFp_leg2 +theme_void() + theme(text = element_text(size=30))


tiff(filename="TEF_mergeF_tree_legend2",height=600,width=800)
print(TEFp_leg2)
dev.off()

#plot&print without legend or labels
TEFp <- plot_tree(TEF_mergeF$tree,type="fan",rank="phylum",taxonomy=TEF_mergeF$taxonomy, size = 1,legend_cutoff =0,main = "TEF")
TEFp <- TEFp +theme_void() + theme(legend.position = "none")


tiff(filename="TEF_mergeF_tree_phylum.tiff",height=650,width=600)
print(TEFp)
dev.off()

#plot&print with legend on right
TEFp_leg2 <- plot_tree(TEF_mergeF$tree,type="fan",rank="phylum",taxonomy=TEF_mergeF$taxonomy, size = 4,legend_cutoff =50,main = "TEF")
TEFp_leg2 <- TEFp_leg2 +theme_void() + theme(text = element_text(size=40))


tiff(filename="TEF_mergeF_tree_phylum_legend2",height=600,width=800)
print(TEFp_leg2)
dev.off()

#plot&print without legend
TEFp <- plot_tree(TEF_mergeF$tree,type="unrooted",rank="phylum",taxonomy=TEF_mergeF$taxonomy, size = 2.5,legend_cutoff =0,main = "TEF")
TEFp <- TEFp +theme_void() + theme(legend.position = "none")


tiff(filename="TEF_mergeF_tree_unrooted_phylum.tiff",height=800,width=800)
print(TEFp)
dev.off()

#plot&print fan tree without legend or labels
TEFp2 <- plot_tree(TEF_mergeF$tree,type="fan",rank="class",taxonomy=TEF_mergeF$taxonomy, size = 1,legend_cutoff =1,main = "TEF") 
TEFp2 <- TEFp2 +theme_void() + theme(legend.position = "none")


tiff(filename="TEF_mergeF_fan_tree",height=650,width=600)
print(TEFp2)
dev.off()

#plot&print phylogram without legend
TEFp3 <- plot_tree(TEF_mergeF$tree,type="phylogram",rank="class",taxonomy=TEF_mergeF$taxonomy, size = 1,legend_cutoff =50,main = "TEF")
TEFp3 <- TEFp3 +theme_void() + theme(legend.position = "none")


tiff(filename="TEF_mergeF_phylo_tree",height=200,width=350)
print(TEFp3)
dev.off()

#### PLOT LENGTHS OF TEF & ITS TOGETHER ####

#create a column with target & merge lengths, then plot

ITS_lengths <- as.data.frame(filtered_lengths_ITS)
colnames(ITS_lengths) <- "Amplicon_length"
ITS_lengths[,2] <- "ITS1"
colnames(ITS_lengths)[2] <- "Target"

TEF_lengths <- as.data.frame(filtered_lengths_TEF)
colnames(TEF_lengths) <- "Amplicon_length"
TEF_lengths[,2] <- "TEF"
colnames(TEF_lengths)[2] <- "Target"

target_lengths <- rbind(TEF_lengths,ITS_lengths)

library(plyr)
mu <- ddply(target_lengths, "Target", summarise, grp.median=median(Amplicon_length))

tiff(filename="Target_lengths_filtered",height=1000,width=750)
ggplot(target_lengths,aes(Amplicon_length,color=Target,fill=Target)) +
  #facet_wrap(~Target) +
  geom_histogram(binwidth=5, alpha=0.5, position="identity") +
  geom_vline(data=mu, aes(xintercept=grp.median, color=Target),
               linetype="dashed") +
  xlab("Amplicon length (bp)") +
  ylab("Count") +
  theme_light() +
  #theme_linedraw() +
  theme(axis.line = element_line(size = 1, colour = "grey20"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=25),)

dev.off()
