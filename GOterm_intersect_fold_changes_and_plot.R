##########################################################################################
##### Goal here is to obtain graphs of GO & Pathway terms that are common to
##### LSS vs ST & LSS vs OSS

##########################################################################################
##### Source & load libraries
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("DESeq2","pasilla","DESeq","pathview"))
# install.packages(c("rJava","ReporteRs","ReporteRsjars","ggplot2","rtable","xtable",
# "VennDiagram","taRifx","devtools","dplyr","data.table"))
# devtools::install_github('tomwenseleers/export',local=F)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R") 
zzz<-c("pheatmap","grid","gplots","ggplot2","export","devtools","DESeq2","pasilla",
       "Biobase","EBSeq","plyr","dplyr","data.table","genefilter","FactoMineR","VennDiagram",
       "DOSE","ReactomePA","org.Hs.eg.db","clusterProfiler","pathview")
lapply(zzz, require, character.only=T)

#########################################################################################
##### Create a function that will plot Terms w/ log2FCs
James_GO_plot <- function(df) {
  title <- deparse(substitute(df)) #substitute is weird (must use it before altering df)
  colnames(df) <- c("V1","LSS_v_ST","LSS_v_OSS")
  df <- melt(df, id.vars=c("V1")) #creates 2 rows from 1 which contains 2 scores
  colnames(df) <- c("Gene","Condition","log2FC")
  df <- arrange(df, desc(log2FC))
  ## set the Gene as a factor, so that ggplot won't plot genenames alphabetically
  df$Gene <- factor(df$Gene, levels = df$Gene)
  ## make the plot & export to PPT
  gg <- ggplot() + geom_bar(aes(y= log2FC, x= Gene, fill= Condition), data= df, stat="identity",
                            position="dodge") + coord_flip() + ggtitle(title) +
    theme(plot.title= element_text(size= 14, face= "bold"),
          axis.text= element_text(size= 14), legend.text= element_text(size= 14),
          legend.title= element_text(size= 14), axis.title= element_text(size= 14))
  print(gg)
  graph2ppt(file="LSS_vs_ST.GO_and_Pathways.pptx", width=4, height=7, append=T)
  dev.off()
}

##### Read in interesting DEGs here
LSS_vs_ST.df <- read.csv(file='../../LSS_vs_ST_x_LSS_vs_OSS.DEGs.csv', quote="\"", header=T)
LSS_vs_ST.GO <- read.delim(file='LSS_vs_ST.ReviGO.GO_bp', quote="\"", header=T)

##### Read interesting GO terms here & merge w/ FC/padj info
### Not that these terms need to be manually ordered as the top 10 in the GO_BP file
blood_vessel_morphogenesis        <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[1,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
extracellular_matrix_organization <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[2,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
regulation_of_cell_migration      <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[3,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)] 
locomotion                        <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[4,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
response_to_wounding              <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[5,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
cell_adhesion                     <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[6,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
wound_healing                     <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[7,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
regulation_of_cell_proliferation  <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[8,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
cell_junction_organization        <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[9,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]
endothelial_cell_differentiation  <- merge(as.data.table(strsplit(as.vector(LSS_vs_ST.GO[10,2]),",")), LSS_vs_ST.df, by.x='V1', by.y='Symbol')[,c(1,3,4)]

##### Run the function on the selected go terms
James_GO_plot(blood_vessel_morphogenesis)
James_GO_plot(extracellular_matrix_organization)
James_GO_plot(regulation_of_cell_migration)
James_GO_plot(locomotion)
James_GO_plot(response_to_wounding)
James_GO_plot(cell_adhesion)
James_GO_plot(wound_healing)
James_GO_plot(regulation_of_cell_proliferation)
James_GO_plot(cell_junction_organization)
James_GO_plot(endothelial_cell_differentiation)

#########################################################################################
##### Create plot of Terms w/ log2FCs for an individual GO TERM
# df <- extracellular_matrix_organization
# title <- deparse(substitute(df)) #substitute is weird (must use it before altering df)
# colnames(df) <- c("V1","LSS_v_ST","LSS_v_OSS")
# df <- melt(df, id.vars=c("V1")) #creates 2 rows from 1 which contains 2 scores
# colnames(df) <- c("Gene","Condition","log2FC")
# df <- arrange(df, desc(log2FC))
# ## set the Gene as a factor, so that ggplot won't plot genenames alphabetically
# df$Gene <- factor(df$Gene, levels = df$Gene)
# ## make the plot & export to PPT
# ggplot() + geom_bar(aes(y= log2FC, x= Gene, fill= Condition), data= df, stat="identity",
#                     position="dodge") + coord_flip() + ggtitle(title) +
#   theme(plot.title= element_text(size= 14, face= "bold"),
#         axis.text= element_text(size= 14), legend.text= element_text(size= 14),
#         legend.title= element_text(size= 14), axis.title= element_text(size= 14))
# graph2ppt(file="LSS_vs_ST_x_LSS_vs_OSS.GO_Pathways.ppt", width=10, height=7, append=T)
