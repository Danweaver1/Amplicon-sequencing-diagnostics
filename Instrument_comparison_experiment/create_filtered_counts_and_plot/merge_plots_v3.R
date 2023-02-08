#### plot merged data ####

file_list=(list.files(pattern="*_merge.csv"))
library(tibble)
library(stringr)
library(ggplot2)
#make colour table assigning each species a specific colour
colour_table <- tibble(
  species = c("Aspergillus_penicillioides","Aspergillus_clavatus","Aspergillus_fischeri","Aspergillus_flavus","Aspergillus_fumigatiaffinis","Aspergillus_fumigatus","Aspergillus_lentulus","Aspergillus_niger","Aspergillus_novofumigatus","Aspergillus_terreus","Candida_albicans","Candida_auris","Candida_glabrata","Candida_inconspicua","Candida_krusei","Candida_parapsilosis","Candida_rugosa","Cryptococcus_neoformans","Fusarium_oxysporum","Homo_sapiens","Lomentospora_prolificans","Penicillium_chrysogenum","Penicillium_rubens","Pichia_kudriavzevii","PjTEF","Rhizopus_arrhizus","Rhizopus_delemar","Rhizopus_oryzae","Saccharomyces_cerevisiae","Aspergillus_novofumigatus","Malassezia_restricta","Meyerozyma_guilliermondii","Aspergillus_nidulans","Aureobasidium_pullulans","Candida_zeylanoides"),
  Colour = c("antiquewhite3", "aquamarine3","azure3","bisque1","blue4","blueviolet","brown4","cadetblue1","chartreuse1","chartreuse4","coral","lavender","darkmagenta","gold2","deeppink1","plum4","grey0","slategray3","tan1","maroon","cyan3","thistle1","yellow","seashell2","blue","honeydew2","tan4","olivedrab4","grey38","#70CA73","#E9E7AE","#E5E5E5","#425865","#599EAC","#C1A3D0"))


colour_table$species_formatted <- stringr::str_replace(colour_table$species, '_', ' ') #create a column for species names without underscores

write.csv(colour_table,"colour_table.csv", row.names=F)  

#check colours
#x11()
#pie(rep(1,length(colour_table$Colour)), col=colour_table$Colour, labels = colour_table$species)

#plot proportional counts
for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,"_merge.csv","")
  plot_name = paste(samplename, ".tiff", sep="")
  MSC=as.data.frame(read.csv(filename, header=FALSE, sep=",", col.names=c("Counts","Species","Sample","Proportional_counts")))
  tiff(plot_name,width = 500, height = 1000)
    print(ggplot(MSC, aes(x=str_replace(Sample,paste(samplename,".", sep=""),""),y=Proportional_counts, fill=Species))
          + scale_fill_manual(breaks = colour_table$species, values= colour_table$Colour, labels=colour_table$species_formatted) #set colours using colour table, set labels using formatted species names
        + geom_bar(stat="identity", width =0.8)
        + scale_y_continuous(limits = c(0, NA),
                               expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
        + ggtitle(samplename)
        + xlab("Sample")
        + ylab("Proportional counts")
        + theme(panel.background = element_rect(fill = "white"), #change background to white
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
                axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
                plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",
                legend.position="bottom")
                + guides(fill=guide_legend(ncol=1))) #sets legend items to one column
  dev.off()
}




  
#plot raw counts

for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,"_merge.csv","")
  plot_name = paste(samplename,"raw", ".tiff", sep="")
  MSC=as.data.frame(read.csv(filename, header=FALSE, sep=",", col.names=c("Counts","Species","Sample","Proportional_counts")))
  tiff(plot_name,width = 500, height = 1000)
  print(ggplot(MSC, aes(x=str_replace(Sample,paste(samplename,".", sep=""),""),y=Counts, fill=Species))
        + scale_fill_manual(breaks = colour_table$species, values= colour_table$Colour, labels=colour_table$species_formatted) #set colours using colour table, set labels using formatted species names
        + geom_bar(stat="identity", width =0.8)
        + scale_y_continuous(limits = c(0, NA),
                             expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
        + ggtitle(samplename)
        + xlab("Sample")
        + ylab("Counts")
        + theme(panel.background = element_rect(fill = "white"), #change background to white
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
                axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
                plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",
                legend.position="bottom")
        + guides(fill=guide_legend(ncol=1))) #sets legend items to one column
  dev.off()
}
 

    