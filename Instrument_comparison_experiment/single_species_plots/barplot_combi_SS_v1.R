library(ggplot2)
library("plyr")
theme_set(theme_classic())
#colnames(otu_mat)  %>% str_remove("unique_species") %>% str_replace("\\.1"," 100 copies") %>% str_replace("\\.2"," 50 copies") %>% str_replace("\\.3"," 20 copies") %>% str_replace("\\.4"," 10 copies") %>% str_replace("\\.5"," 2 copies")  %>% str_remove("_LOD[0-9]") %>% str_replace("TEF[0-9]*_","TEF ")  %>% str_replace("ITS[0-9]*_","ITS1 ")

#### plot merged data barplots ####

file_list=(list.files(pattern="*_merge.csv"))
library(tibble)
library(stringr)


#make colour table assigning each species a specific colour
#make colour table assigning each species a specific colour
colour_table <- tibble(
  species = c("Aspergillus_penicillioides","Aspergillus_penicilloides","Aspergillus_clavatus","Aspergillus_fischeri","Aspergillus_flavus","Aspergillus_fumigatiaffinis","Aspergillus_fumigatus","Aspergillus_niger","Aspergillus_niger-or-welwitschiae-or-lacticoffeatus-or-awamori","Aspergillus_novofumigatus-or-udagawae-or-lentulus","Aspergillus_terreus","Candida_albicans","Candida_auris","Candida_glabrata","PjTEF","Candida_inconspicua","Candida_krusei","Candida_rugosa","Cryptococcus_neoformans","Fusarium_oxysporum","Homo_sapiens","Candida_parapsilosis","Penicillium_chrysogenum","Penicillium_rubens","Pichia_kudriavzevii","Lomentospora_prolificans","Rhizopus_arrhizus","Rhizopus_delemar","Rhizopus_oryzae","Saccharomyces_cerevisiae","Malassezia_restricta","Meyerozyma_guilliermondii","Aureobasidium_pullulans","Aspergillus_lentulus"),
  Colour = c("antiquewhite3","antiquewhite3", "brown4","azure2","bisque1","blue4","darkorchid4","cadetblue","cadetblue1","chartreuse1","chartreuse4","darkorange2","lavender","darkmagenta","gold2","deeppink1","plum4","grey0","azure3","tan1","maroon","deepskyblue","yellow","thistle1","seashell2","seagreen4","honeydew2","tan4","olivedrab4","grey38","paleturquoise3","#FF5151","#B95C2E","#FDE56B"))






colour_table$species_formatted <- stringr::str_replace(colour_table$species, '_', ' ') #create a column for species names without underscores

write.csv(colour_table,"colour_table.csv", row.names=F)  


#plot proportional counts
for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,"_merge.csv","")
  plot_name = paste(samplename, ".tiff", sep="")
  MSC=as.data.frame(read.csv(filename, header=FALSE, sep=",", col.names=c("Counts","Species","Sample","Proportional_counts")))
  MSC$Sample <- MSC$Sample %>% str_replace('S0[0-9].', '') %>% str_replace('3_ITS1', 'ITS1') %>% str_replace('3R_ITS1', 'ITS1') %>% str_replace("_MiSeq","") %>% str_replace("3_TEF","TEF BC 1") %>% str_replace("1_TEF","TEF BC 3") %>% str_replace("2_TEF","TEF BC 2")
  tiff(plot_name,width = 700, height = 600)
  print(ggplot(MSC, aes(x=str_replace(Sample,paste(samplename,".", sep=""),""),y=Proportional_counts, fill=Species)) +
      #set colours using colour table, set labels using formatted species names 
      geom_bar(stat="identity", position="fill", width =0.8) +
      scale_fill_manual(breaks = colour_table$species, values= colour_table$Colour, labels=colour_table$species_formatted) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) + #makes Y axis start from origin and expands so that 1 is not cut off
      ggtitle(samplename) +
      xlab("Sample") +
      ylab("Proportion") +
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
                    legend.position="right") +
      guides(fill=guide_legend(ncol=1,keywidth = 2)))
  dev.off()#sets legend items to one column
}


#plot raw counts

for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,"_merge.csv","")
  plot_name = paste(samplename, "raw.tiff", sep="")
  MSC=as.data.frame(read.csv(filename, header=FALSE, sep=",", col.names=c("Counts","Species","Sample","Proportional_counts")))
  MSC$Sample <- MSC$Sample %>% str_replace('S0[0-9].', '') %>% str_replace('3_ITS1', 'ITS1') %>% str_replace('3R_ITS1', 'ITS1') %>% str_replace("_MiSeq","") %>% str_replace("3_TEF","TEF BC 1") %>% str_replace("1_TEF","TEF BC 3") %>% str_replace("2_TEF","TEF BC 2")
  tiff(plot_name,width = 1200, height = 800)
  print(ggplot(MSC, aes(x=str_replace(Sample,paste(samplename,".", sep=""),""),y=Counts, fill=Species)) +
          #set colours using colour table, set labels using formatted species names 
          geom_bar(stat="identity", width =0.8) +
          scale_fill_manual(breaks = colour_table$species, values= colour_table$Colour, labels=colour_table$species_formatted) +
          scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) + #makes Y axis start from origin and expands so that 1 is not cut off
          ggtitle(samplename) +
          xlab("Sample") +
          ylab("Counts") +
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
                legend.position="right") +
          guides(fill=guide_legend(ncol=1,keywidth = 2)))
  dev.off()#sets legend items to one column
}




