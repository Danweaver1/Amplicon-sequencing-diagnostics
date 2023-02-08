
df <- read.csv("ID_rate.csv",header=TRUE,row.names=1)
#get working dir
dir <- getwd()
#take basename (ie folder actually in)
basedir <- basename(dir)
#create table with ID rates only
x <- df[,((length(df) - length(df)/3)+1):length(df)]
ID_rates_only <- subset(x, rownames(x) != "PjTEF")
library("pheatmap")
library("stringr")

#plot & send heatmap to file

png("ID_rates_w_vals.png",width = 800, height = 800)
  pheatmap(ID_rates_only, display_numbers = trunc(ID_rates_only, prec=4), color = colorRampPalette(c("firebrick3", "white", "chartreuse4"))(25),fontsize=16, cellwidth=30,main=basedir)
dev.off()




colnames(ID_rates_only) <- colnames(ID_rates_only) %>% str_replace_all("[_]", " ") %>% str_replace("PE ","")  %>% str_replace("merge ","") %>% str_replace("ID rate","") %>% str_replace("IDrate","")  %>% str_replace("ITS","ITS1")
rownames(ID_rates_only) <- rownames(ID_rates_only) %>% str_replace_all("[_]", " ")

newnames <- lapply(
  rownames(ID_rates_only),
  function(x) bquote(italic(.(x))))

png("ID_rates_w_vals_Tmerge_Ipe.png",width = 800, height = 1000)
pheatmap(ID_rates_only, display_numbers = trunc(ID_rates_only, prec=4), number_color = "grey20", cluster_rows = FALSE,
         cluster_cols = FALSE, color = colorRampPalette(c("#de2d26","#f7f7f7","#67a9cf"))(25),fontsize=20, cellwidth=50,cellheight = 40,labels_row = as.expression(newnames))
dev.off()

########plot with annotation for target and sequencer
my_sample_col <- data.frame(Target = rep(c("ITS1", "TEF"), c(3,2)),Sequencer=c(rep("iSeq",2), "MiSeq","iSeq","MiSeq"))
row.names(my_sample_col) <- colnames(ID_rates_only)

my_colour = list(
  Target = c(ITS1 = "#F8766D", TEF = "#00BFC4"),
  Sequencer = c(iSeq = "#ccebc5", MiSeq = "#1b9e77")
)
png("ID_rates_w_vals_annot.png",width = 800, height = 1000)
pheatmap(ID_rates_only, display_numbers = trunc(ID_rates_only, prec=4), number_color = "grey20", cluster_rows = FALSE,
         cluster_cols = FALSE, annotation_col = my_sample_col, annotation_colors = my_colour, color = colorRampPalette(c("#de2d26","#f7f7f7","#67a9cf"))(25),fontsize=15, cellwidth=30,cellheight = 25,labels_row = as.expression(newnames))
dev.off()


#reorder species names to fit with FC plot

new_ID_rates_only <- ID_rates_only[ order(row.names(ID_rates_only),decreasing = TRUE), ]

newnames2 <- lapply(
  rownames(new_ID_rates_only),
  function(x) bquote(italic(.(x))))
png("ID_rates_w_vals_Tmerge_Ipe_reorder.png",width = 800, height = 1000)
pheatmap(new_ID_rates_only, display_numbers = trunc(new_ID_rates_only, prec=4), number_color = "grey20", cluster_rows = FALSE,
         cluster_cols = FALSE, color = colorRampPalette(c("#de2d26","#f7f7f7","#67a9cf"))(25),fontsize=20, cellwidth=50,cellheight = 40,labels_row = as.expression(newnames2))
dev.off()
#display_numbers = trunc(ID_rates_only, prec=4) - adds the actually value to each square
#color = colorRampPalette(c("firebrick3", "white", "chartreuse4"))(25) - sets colour ramp to have 25 distinct colours and go from red to white to green (green being 100%)
