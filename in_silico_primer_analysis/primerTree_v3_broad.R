## Performs in silico analysis of primer hits using primerTree package ####
##  No restriction is placed on the organisms searched 
## NOTE! Requires an NCBI API key to run

#### setup ####
#install.packages("primerTree")
library(primerTree)
library(ape)
library(ggplot2)
library(dplyr)
library(rstatix)

#### ITS1 degenerate primer set, 1000 aligns, 50 permutations ####
ITS1_broad <- search_primer_pair(name='ITS1','TCCGTAGGTGAACCTGCGG', 'GGCYRCGTTCTTCATCGAYGC',
                                 total_primer_specificity_mismatch=3,num_aligns=1000,num_permutations=50,
                                 api_key="")

tiff(filename="ITS1_broad",height=1000,width=1000)
print(plot(ITS1_broad))
dev.off()

ITS1_broad_lengths <- seq_lengths(ITS1_broad,summarize=TRUE)

tiff(filename="ITS1_broad_lengths",height=1000,width=1000)
print(barplot(ITS1_broad_lengths,main = "Frequency of sequence lengths for ITS1 primers",xlab="Amplicon length (in bp)",ylab=("Frequency")))
dev.off()

ITS1_broad_Blast <- ITS1_broad$BLAST_result
write.csv(ITS1_broad_Blast,file="ITS1_broad_blast_results")

length(ITS1_broad_Blast$gi)
length(unique(ITS1_broad_Blast$gi))

median(ITS1_broad_Blast$product_length)
mean(ITS1_broad_Blast$product_length)
range(ITS1_broad_Blast$product_length)

ITS1_broad_summ <- get_summary_stats(ITS1_broad_Blast)
write.csv(ITS1_broad_summ,file="ITS1_broad_results_summary")


write.csv(ITS1_broad$taxonomy,file="ITS1_broad_taxonomy")
write.tree(ITS1_broad$tree,file="ITS1_broad_tree")

##plot trees ##

#plot&print with legend on bottom
ITSp_leg <- plot_tree(ITS1_broad$tree,type="radial",rank="kingdom",taxonomy=ITS1_broad$taxonomy, size = 5,legend_cutoff =50,main = "ITS1")
ITSp_leg <- ITSp_leg +theme_void() + theme(text = element_text(size=30),legend.position = "bottom")

tiff(filename="ITS1_broad_tree_legend1",height=1000,width=1600)
print(ITSp_leg)
dev.off()



#plot&print KINGDOM unrooted tree without legend or labels
ITSp2 <- plot_tree(ITS1_broad$tree,type="fan",rank="kingdom",taxonomy=ITS1_broad$taxonomy, size = 2.5,legend_cutoff =1,main = "ITS1") 
ITSp2 <- ITSp2 +theme_void() + theme(legend.position = "none")


tiff(filename="ITS1_broad_kingdom_fan_tree",height=800,width=800)
print(ITSp2)
dev.off()


####TEF merge  ####

TEF_merge_broad <- search_primer_pair(name='TEF_merge_broad','YGGYGARTTCGARGCYGGTAT',
                                      'ADCATRTTRTCDCCRTKGMADCC',total_primer_specificity_mismatch=3,
                                      num_aligns=1000,num_permutations=50,
                                      api_key="")


tiff(filename="TEF_merge_broad",height=1000,width=1000)
print(plot(TEF_merge_broad))
dev.off()

TEF_merge_broad_lengths <- seq_lengths(TEF_merge_broad,summarize=TRUE)

tiff(filename="TEF_merge_broad_lengths",height=1000,width=1000)
print(barplot(TEF_merge_broad_lengths,main = "Frequency of sequence lengths for TEF primers",xlab="Amplicon length (in bp)",ylab=("Frequency")))
dev.off()

TEF_merge_broad_Blast <- TEF_merge_broad$BLAST_result
write.csv(TEF_merge_broad_Blast,file="TEF_merge_broad_blast_results")

length(TEF_merge_broad_Blast$gi)
length(unique(TEF_merge_broad_Blast$gi))

median(TEF_merge_broad_Blast$product_length)
mean(TEF_merge_broad_Blast$product_length)
range(TEF_merge_broad_Blast$product_length)



TEF_merge_broad_summ <- get_summary_stats(TEF_merge_broad_Blast)
write.csv(TEF_merge_broad_summ,file="TEF_merge_broad_results_summary")


write.csv(TEF_merge_broad$taxonomy,file="TEF_merge_broad_taxonomy")
write.tree(TEF_merge_broad$tree,file="TEF_merge_broad_tree")


##plot trees ##

#plot&print with legend on bottom
TEFp_leg <- plot_tree(TEF_merge_broad$tree,type="radial",rank="kingdom",taxonomy=TEF_merge_broad$taxonomy, size = 5,legend_cutoff =50,main = "TEF")
TEFp_leg <- TEFp_leg +theme_void() + theme(text = element_text(size=30),legend.position = "bottom")

tiff(filename="TEF_merge_broad_tree_legend1",height=1000,width=1600)
print(TEFp_leg)
dev.off()



#plot&print KINGDOM unrooted tree without legend or labels
TEFp2 <- plot_tree(TEF_merge_broad$tree,type="fan",rank="kingdom",taxonomy=TEF_merge_broad$taxonomy, size = 2.5,legend_cutoff =1,main = "TEF") 
TEFp2 <- TEFp2 +theme_void() + theme(legend.position = "none")


tiff(filename="TEF_merge_broad_kingdom_fan_tree",height=800,width=800)
print(TEFp2)
dev.off()
