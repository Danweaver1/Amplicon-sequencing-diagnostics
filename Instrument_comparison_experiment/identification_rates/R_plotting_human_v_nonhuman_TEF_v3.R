
setwd("./TEF_human_10X")
library("stringr")
df <- read.csv("ID_rate.csv",header=TRUE,row.names=1)

#merge iseq and miseq data 

df[,"merge_TEF_ID"] <- df[,"merge_TEF_iSeq_ID"] + df[,"merge_TEF_MiSeq_ID"]
df[,"merge_TEF_Missing"] <- df[,"merge_TEF_iSeq_Missing"] + df[,"merge_TEF_MiSeq_Missing"]
#recalculate %ID

df[,"merge_TEF_ID_rate"] <- (df[,"merge_TEF_ID"]/(df[,"merge_TEF_ID"] + df[,"merge_TEF_Missing"]))*100
df[,"merge_TEF_ID_rate"]


#reduce and add + to indicate human background to colnames
library(stringr)
colnames(df) <- colnames(df) %>% str_replace_all("[_]", " ") %>% str_replace("PE ","")  %>% str_replace("merge ","") %>% str_replace("ID rate","+") %>% str_replace("IDrate","+")



#read in ID rates for equivalent samples without human background
setwd("../TEF_no_human")
df1 <- read.csv("ID_rate.csv",header=TRUE,row.names=1)

#merge iseq and miseq data 

df1[,"merge_TEF_ID"] <- df1[,"merge_TEF_iSeq_ID"] + df1[,"merge_TEF_MiSeq_ID"]
df1[,"merge_TEF_Missing"] <- df1[,"merge_TEF_iSeq_Missing"] + df1[,"merge_TEF_MiSeq_Missing"]
#recalculate %ID
df1[,"merge_TEF_ID_rate"] <- (df1[,"merge_TEF_ID"]/(df1[,"merge_TEF_ID"] + df1[,"merge_TEF_Missing"]))*100

colnames(df1) <- colnames(df1) %>% str_replace_all("[_]", " ") %>% str_replace("PE ","")  %>% str_replace("merge ","") %>% str_replace("ID rate","-") %>% str_replace("IDrate","-")


full_tbl <- cbind(df[,8:9],df1[,8:9])
#create table with combined ID rates only
x <- full_tbl[,c(2,4)]

ID_rates_only <- subset(x, rownames(x) != "PjTEF")

setwd('..')
library("pheatmap")


#improve plot text and plot only PE ITS and merge TEF



rownames(ID_rates_only) <- rownames(ID_rates_only) %>% str_replace_all("[_]", " ")

#reorder sample cols
ID_rates_only <- ID_rates_only[ ,order(colnames(ID_rates_only)) ]
newnames <- lapply(
  rownames(ID_rates_only),
  function(x) bquote(italic(.(x))))

png("ID_rates_w_vals.png",width = 800, height = 1000)
pheatmap(ID_rates_only, display_numbers = trunc(ID_rates_only, prec=4), number_color = "grey20", cluster_rows = FALSE,
         cluster_cols = FALSE, color = colorRampPalette(c("#de2d26","#f7f7f7","#67a9cf"))(25),fontsize=20, cellwidth=50,cellheight = 40,labels_row = as.expression(newnames))
dev.off()



########plot with annotation for target and sequencer
my_sample_col <- data.frame("Human_background"=c("No", "Yes"))
row.names(my_sample_col) <- colnames(ID_rates_only)

my_colour = list(
  "Human_background" = c(No = "#999999", Yes = "#E69F00"))


png("ID_rates_w_vals_annot.png",width = 800, height = 1000)
pheatmap(ID_rates_only, display_numbers = trunc(ID_rates_only, prec=4), number_color = "grey20", cluster_rows = FALSE,
         cluster_cols = FALSE, annotation_col = my_sample_col, annotation_colors = my_colour, color = colorRampPalette(c("#de2d26","#f7f7f7","#67a9cf"))(25),fontsize=20, cellwidth=40,cellheight = 25,labels_row = as.expression(newnames))
dev.off()


#reorder species names to fit with FC plot


new_ID_rates_only <- ID_rates_only[ order(row.names(ID_rates_only),decreasing = TRUE), ]

newnames2 <- lapply(
  rownames(new_ID_rates_only),
  function(x) bquote(italic(.(x))))
png("ID_rates_w_vals_reorder.png",width = 800, height = 1000)
pheatmap(new_ID_rates_only, display_numbers = trunc(new_ID_rates_only, prec=4), number_color = "grey20", cluster_rows = FALSE,
         cluster_cols = FALSE, color = colorRampPalette(c("#de2d26","#f7f7f7","#67a9cf"))(25),fontsize=20, cellwidth=50,cellheight = 40,labels_row = as.expression(newnames2))
dev.off()
#display_numbers = trunc(ID_rates_only, prec=4) - adds the actually value to each square
#color = colorRampPalette(c("firebrick3", "white", "chartreuse4"))(25) - sets colour ramp to have 25 distinct colours and go from red to white to green (green being 100%)