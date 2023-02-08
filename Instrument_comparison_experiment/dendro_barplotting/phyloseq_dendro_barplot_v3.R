### Phyloseq 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("phyloseq")

library("phyloseq"); packageVersion("phyloseq")
library("plyr")
library("ggplot2"); packageVersion("ggplot2")
library("stringr")
theme_set(theme_bw())
theme_set(theme_classic())
samplename="MD1"
#or: theme_set(theme_classic()) or whatever
##Run this script on output of join TEF v1 shell script - feed in unique species counts with filename format: targetsamplenumber_instrument_mix_unique_species
# Read files in

otu_mat <- as.matrix(read.table("out", sep ="\t", header=TRUE))
colnames(otu_mat)<- colnames(otu_mat)  %>% str_replace("rep2","B") %>% str_remove("[0-9][0-9]") %>%  str_remove("M[A-Z][0-9]") %>% str_remove("unique_species") %>% str_replace("_", " ") %>% str_remove_all("_") %>% str_replace(".[0-9] Mock","Mock") %>% str_replace(".1"," rep1") %>% str_replace(".2"," rep2") %>% str_replace("ITS","ITS1")

tax_mat <- as.matrix(read.table("fungi_tax_TEF_v4.tsv", sep ="\t", header=TRUE))

# Define row names from otu column "X"
row.names(otu_mat) <- otu_mat[,1]

# Remove the comun from the matrix
otu_mat <- otu_mat[,-1]

# Save the data as numeric
class(otu_mat) <- "numeric"


# Repeat for tax
row.names(tax_mat) <- tax_mat[,1]
colnames(tax_mat)

# Transform into phyloseq objects
OTU3 = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX3 = tax_table(tax_mat)
physeq1 = phyloseq(OTU3, TAX3) 

##add sample info
#target or mock?
samp_types <- sample_names(physeq1) %>% str_remove(" MiSeq")  %>% str_remove(" iSeq")  %>% str_remove(" rep[0-9]") %>% str_remove("B")

class(samp_types)
#create sample data 
SAMP = sample_data(data.frame(
 SampleType =samp_types,
  row.names=sample_names(physeq1),
  stringsAsFactors=TRUE
))

#merge sample data with phyloseq object
physeq = merge_phyloseq(physeq1, SAMP)


#normalise counts
total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeqn = transform_sample_counts(physeq, standf)
physeq_abund <- filter_taxa(physeqn, function(x) sum(x > total*0.002) > 0, TRUE)

###check what remains after normalising
#ntaxa(physeq_abund)
#taxa_names(physeq_abund)
#nsamples(physeq_abund)
#sample_names(physeq_abund)
#sample_variables(physeq_abund)
#otu_table(physeq_abund)

###use physeq distance function to calclate distance using bray
##see all possible methods (using bray here)
#dist_methods <- unlist(distanceMethodList)
#print(dist_methods)
##calc distance
d <-distance(physeq_abund,method="bray",type="samples")

write.csv(as.matrix(d),file=paste(samplename,"_distances.csv",sep=""))
#cluster distances and create dendrogram
dend <- as.dendrogram(hclust(d))
order.dendrogram(dend)
str(dend)
#change dendro order so that mock is in middle

#dend.reorder <- reorder(dend,c(6,7,8,9,5,4:1))
dend.reorder <- reorder(dend,c(11:9,7,8,6:1))
str(dend.reorder)



#convert normalised otu table to lnog format for plotting
library(tidyr)
#extract otu table from physeq abunbance object and convert to df
norm_otu <- as.data.frame(otu_table(physeq_abund))
#add species column using rownames
norm_otu[(length(norm_otu)+1)] <- row.names(norm_otu)
colnames(norm_otu)[length(norm_otu)] <- "Species"
#make long
data_long <- gather(norm_otu, Sample, Normalised_counts, -c(Species))


library(ggplot2)
library(ggdendro)
library(dendextend)
library(cowplot)

##this line is what orders barplot same as dendro (was already in v1 script)
data_long$Sample <- factor(data_long$Sample, levels = labels(dend.reorder))

#First, plot dendrogram and barplots together to check order of samples in dendro
##need to reorder samples in barplot to match dendro
p2 <- ggplot(data_long, aes(fill=Species, y=Normalised_counts, x=Sample)) + 
  geom_bar( stat="identity", position="fill") +
  theme(legend.position = "none", axis.title.y=element_blank())


p1 <- ggdendrogram(dend.reorder,rotate=F,labels=TRUE,leaf_labels =TRUE)

#+ coord_flip()

#plot_grid(p1, p2, align = "h")
plot_grid(p1, p2, ncol=1)
#send to file
dend_bar_plotname <- paste(samplename, "norm_abund_dendro.tiff",sep="")
tiff(dend_bar_plotname,width = 500, height = 600)
print(plot_grid(p1, p2, ncol=1))
dev.off()


#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

##can extract data to get the values for segments
dend_data <- dendro_data(dend.reorder, type = "rectangle")
names(dend_data)
dend_data$segments
dend_data$labels

segs <- dend_data$segments
#find range of values so can alter axis limits
range(segs$y)
range(segs$yend)
# Plot line segments and add labels (plot dendrogram alone, for manaully adding to barplot for figure)
#x11()
ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = " "),
            hjust = 1, angle = 90, size = 3)+ theme_dendro() +
  ylim(-0.1, 0.4)

dend_plotname <- paste(samplename,"dendro", ".tiff", sep="")
tiff(dend_plotname,width = 800, height = 600)
print(ggplot(dend_data$segments) + 
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
        geom_text(data = dend_data$labels, aes(x, y, label = " "),
                  hjust = 1, angle = 90, size = 3)+ theme_dendro() +
        ylim(-0.1, 3))
dev.off()


#### plot merged data barplots ####

library(tibble)
library(stringr)
library(ggplot2)
#make colour table assigning each species a specific colour
colour_table <- tibble(
  species = c("Aspergillus_penicillioides","Aspergillus_penicilloides","Aspergillus_clavatus","Aspergillus_fischeri","Aspergillus_flavus","Aspergillus_fumigatiaffinis","Aspergillus_fumigatus","Aspergillus_niger","Aspergillus_niger-or-welwitschiae-or-lacticoffeatus-or-awamori","Aspergillus_novofumigatus-or-udagawae-or-lentulus","Aspergillus_terreus","Candida_albicans","Candida_auris","Candida_glabrata","PjTEF","Candida_inconspicua","Candida_krusei","Candida_rugosa","Cryptococcus_neoformans","Fusarium_oxysporum","Homo_sapiens","Candida_parapsilosis","Penicillium_chrysogenum","Penicillium_rubens","Pichia_kudriavzevii","Lomentospora_prolificans","Rhizopus_arrhizus","Rhizopus_delemar","Rhizopus_oryzae","Saccharomyces_cerevisiae"),
  Colour = c("antiquewhite3","antiquewhite3", "brown4","azure3","bisque1","blue4","darkorchid4","cadetblue","cadetblue1","chartreuse1","chartreuse4","darkorange2","lavender","darkmagenta","gold2","deeppink1","plum4","grey0","azure3","tan1","maroon","deepskyblue","thistle1","yellow","seashell2","seagreen4","honeydew2","tan4","olivedrab4","grey38"))


colour_table$species_formatted <- stringr::str_replace(colour_table$species, '_', ' ') #create a column for species names without underscores

write.csv(colour_table,"colour_table.csv", row.names=F)  

#plot Normalised counts


plot_name = paste(samplename, ".tiff", sep="")

tiff(plot_name,width = 500, height = 900)
print(ggplot(data_long, aes(x=Sample,y=Normalised_counts, fill=Species)) +
    #set colours using colour table, set labels using formatted species names 
    geom_bar(stat="identity", position="fill", width =0.8) +
    scale_fill_manual(breaks = colour_table$species, values= colour_table$Colour, labels=colour_table$species_formatted) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) + #makes Y axis start from origin and expands so that 1 is not cut off
    ggtitle(samplename) +
    xlab("Sample") +
    ylab("Normalised counts") +
    theme(panel.background = element_rect(fill = "white"), #change background to white
                  axis.line = element_line(size = 1, colour = "black"),
                  axis.ticks = element_line(size = 2),
                  axis.text.x = element_text(color = "grey20", size = 20, angle =65, hjust = 1, vjust = 1, face = "plain"),
                  axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                  axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
                  axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
                  plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
                  legend.key.size = unit(3, 'lines'),
                  legend.title=element_text(size=27),
                  legend.text=element_text(size=25, face="italic"),
                  legend.direction="vertical",
                  legend.position="bottom") +
    guides(fill=guide_legend(ncol=1,keywidth = 2)))
dev.off()#sets legend items to one column





