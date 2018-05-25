################ Integrate RNAseq DEGs with ATAC-seq DA peaks ####################
##################################################################################
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R") 
# source("https://bioconductor.org/biocLite.R")
# biocLite(pkgs = c("pheatmap","gplots","ggplot2","devtools","DESeq2","pasilla",
# "Biobase","EBSeq","dplyr","data.table","genefilter","FactoMineR","VennDiagram","DOSE","ReactomePA",
# "org.Hs.eg.db","clusterProfiler","pathview","ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene",
# "GO.db", "GeneOverlap", "GOSemSim"))
zzz<-c("pheatmap","grid","gplots","ggplot2","export","devtools","DESeq2","pasilla","Biobase","EBSeq",
       "dplyr","data.table", "genefilter","FactoMineR","VennDiagram","DOSE","ReactomePA","org.Hs.eg.db",
       "clusterProfiler","pathview","ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","GO.db",
       "GeneOverlap","GOSemSim")
lapply(zzz, require, character.only= T)
# devtools::install_github('tomwenseleers/export',local=F) # to re-install export, if necessary
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

##################################################################################
##### import Lam_vs_Static DA peaks & subsetted up/down regulated
# DA_sites_w_FC <- read.delim("Lam_vs_Static.edgeR_blocked.DA_sites")
LSS_enriched_DA.proximal_genes.df <- read.delim("LSS_enriched_DA.proximal_genes.txt", quote="\"'")
ST_enriched_DA.proximal_genes.df <- read.delim("ST_enriched_DA.proximal_genes.txt", quote="\"'")
LSS_enriched_DA.proximal_gene.EntrezIDs <- as.character(LSS_enriched_DA.proximal_genes.df$geneId)
ST_enriched_DA.proximal_genes.EntrezIDs <- as.character(ST_enriched_DA.proximal_genes.df$geneId)

##### Read in the name conversion file
name_convert <- read.delim('kgXref_NO_EXCEL.txt', quote="\"'")

##### merge peaks with gene Symbols & bind to make combined DA peakset
up_peaks_symbol<- merge(LSS_enriched_DA.proximal_genes.df, name_convert, by.x= "transcriptId", by.y= "transcriptId")
down_peaks_symbol<- merge(ST_enriched_DA.proximal_genes.df, name_convert, by.x= "transcriptId", by.y= "transcriptId")
up_and_down_DApeaks <- as.data.table(rbind(up_peaks_symbol, down_peaks_symbol))

##### merge all DA proximal peaks w/ associated FC values
LSS_vs_ST.DA_w_FCs <- as.data.table(read.delim('1701.LSS_vs_ST.edgeR_blocked.DA_sites', quote="\"'"))
setkey(LSS_vs_ST.DA_w_FCs, seqnames, start, end)
LSS_vs_ST.proximal_DA_w_FC <- foverlaps(up_and_down_DApeaks, LSS_vs_ST.DA_w_FCs, type="any")
write.table(LSS_vs_ST.proximal_DA_w_FC, "1701.LSS_vs_ST.proximal_DA_w_FC.txt", sep="\t", row.names= F)

##### import Lam_vs_Static DEGs
ST_v_LSS.DEGs <- read.delim("1701.ST_v_LSS.FDR10.EBSeq_results")

##### subset LSS enriched DEGs & obtain EntrezIDs
DEGs.LSS.df <- subset(ST_v_LSS.DEGs, log2FC.1701.Static_v_Lam < 0) #subset DEGs up in Lam
DEGs.LSS.Ensembl <- as.vector(DEGs.LSS.df$Ensembl_nodot)
intersected.Ens_Entrez_LSS <- clusterProfiler:: bitr(DEGs.LSS.Ensembl, fromType="ENSEMBL",
                                                     toType="ENTREZID", OrgDb="org.Hs.eg.db")
DEGs.LSS.EntrezID <- as.vector(intersected.Ens_Entrez_LSS[['ENTREZID']])

##### subset ST enriched DEGs & obtain EntrezIDs
DEGs.ST.df <- subset(ST_v_LSS.DEGs, log2FC.1701.Static_v_Lam > 0) #subset DEGs up in Lam
DEGs.ST.Ensembl <- as.vector(DEGs.ST.df$Ensembl_nodot)
intersected.Ens_Entrez_ST <- clusterProfiler:: bitr(DEGs.ST.Ensembl, fromType="ENSEMBL",
                                                     toType="ENTREZID", OrgDb="org.Hs.eg.db")
DEGs.ST.EntrezID <- as.vector(intersected.Ens_Entrez_ST[['ENTREZID']])

#############################################################################################
############### Test/plot correlations between proximal DA peaks & associated gene expression
LvS.prox_DA_FC <- LSS_vs_ST.proximal_DA_w_FC[,c(1:4,9,27,28),with=F] ## Subset desired columns
all_gene_expression <- read.delim('?not sure yet', quote="\"'")

## merge peaks w all_gene_expression & plot both & regress for correlation

#############################################################################################
#####################  investigate DEG and peak overlaps with venn diagrams
## overlap all DEGs with all promoter DAs
venn.plot.DEG_x_DA <- venn.diagram(list(DEGs= as.vector(Static_v_Lam.DEGs[['Symbol']]),
                                            DAs= as.vector(up_and_down_DApeaks[['Symbol']])),
                                       filename= NULL, main= "Laminar vs Static DEGs X DA")
grid.draw(venn.plot.DEG_x_DA)
graph2ppt(file="DEG_vs_DA.ppt", width=2, height=2, append=T)
dev.off()
##
venn.plot.DEGup_x_DAup <- venn.diagram(list(DEG_up= as.vector(up_DEGs[['Symbol']]),
                                        DA_up= as.vector(up_peaks_symbol[['Symbol']])),
                                   filename= NULL, main= "Laminar vs Static DEGs X DA")
grid.draw(venn.plot.DEGup_x_DAup)
graph2ppt(file="DEG_vs_DA.ppt", width=2, height=2, append=T)
dev.off()
#
venn.plot.DEGdown_x_DAdown <- venn.diagram(list(DEG_down= as.vector(down_DEGs[['Symbol']]),
                                            DA_down= as.vector(down_peaks_symbol[['Symbol']])),
                                       filename= NULL, main= "Laminar vs Static DEGs X DA")
grid.draw(venn.plot.DEGdown_x_DAdown)
graph2ppt(file="DEG_vs_DA.ppt", width=2, height=2, append=T)
dev.off()
#
venn.plot.DEGup_x_DAdown <- venn.diagram(list(DEG_up= as.vector(up_DEGs[['Symbol']]),
                                            DA_down= as.vector(down_peaks_symbol[['Symbol']])),
                                       filename= NULL, main= "Laminar vs Static DEGs X DA")
grid.draw(venn.plot.DEGup_x_DAdown)
graph2ppt(file="DEG_vs_DA.ppt", width=2, height=2, append=T)
dev.off()
#
venn.plot.DEGdown_x_DAup <- venn.diagram(list(DEG_down= as.vector(down_DEGs[['Symbol']]),
                                            DA_up= as.vector(up_peaks_symbol[['Symbol']])),
                                       filename= NULL, main= "Laminar vs Static DEGs X DA")
grid.draw(venn.plot.DEGdown_x_DAup)
graph2ppt(file="DEG_vs_DA.ppt", width=2, height=2, append=T)
dev.off()

##### Merge DA_up & DEG_up to get intersected list
upDA_upDEG <- merge(up_DEGs, up_peaks_symbol, by.x= "Symbol", by.y= "Symbol")
write.table(upDA_upDEG, "Lam_v_static.upDA_upDEG.txt", sep="\t", row.names= F)

##### Merge DA_down & DEG_down to get intersected list
downDA_downDEG <- merge(down_DEGs, down_peaks_symbol, by.x= "Symbol", by.y= "Symbol")
write.table(downDA_downDEG, "Lam_v_static.downDA_downDEG.txt", sep="\t", row.names= F)

##########################################################################################
##### GeneOverlap Package for significance testing of overlaps

## DA_up & DEG_up
go.obj.up <- newGeneOverlap(up_peaks_symbol$Symbol, up_DEGs$Symbol, genome.size=57820)
go.obj.up <- testGeneOverlap(go.obj.up)
go.obj.up
## DA_down & DEG_down
go.obj.down <- newGeneOverlap(down_peaks_symbol$Symbol, down_DEGs$Symbol, genome.size=57820)
go.obj.down <- testGeneOverlap(go.obj.down)
go.obj.down

## Get p-values & odds.ratio's for all DEG vs DA comparisons
DA_up <- unique(as.character(up_peaks_symbol$Symbol))
DA_down <- unique(as.character(down_peaks_symbol$Symbol))
DEG_up <- as.character(up_DEGs$Symbol)
DEG_down <- as.character(down_DEGs$Symbol)
Peaks <- list(DA_up=DA_up,DA_down=DA_down)
DEGs <- list(DEG_up=DEG_up,DEG_down=DEG_down)
gom.obj <- newGOM(Peaks, DEGs, genome.size=57820)
## colorkey represents the odds ratios & sig p-values are superimposed on the grids
drawHeatmap(gom.obj, what="odds.ratio", ncolused=4, grid.col="Oranges", note.col="black")
graph2ppt(file="DEG_vs_DA.ppt", width=6, height=6, append=T)
dev.off()

##########################################################################################
##### GoSemSim to compare semantic similarity of GO terms between conditions

GO_genes <- list(LSS_DA=   c(LSS_enriched_DA.proximal_gene.EntrezIDs),
                 ST_DA=    c(ST_enriched_DA.proximal_genes.EntrezIDs),
                 LSS_DEGs= c(DEGs.LSS.EntrezID),
                 ST_DEGs=  c(DEGs.ST.EntrezID))

##### Test for similarity in combined GO terms sets (MF)
LSS_GO_MF_SemSim <- clusterSim(GO_genes[["LSS_DEGs"]], GO_genes[["LSS_DA"]], ont= "MF",
                              organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
ST_GO_MF_SemSim <- clusterSim(GO_genes[["ST_DEGs"]], GO_genes[["ST_DA"]], ont= "MF",
                              organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
DEGs_GO_MF_SemSim <- clusterSim(GO_genes[["LSS_DEGs"]], GO_genes[["ST_DEGs"]], ont= "MF",
                              organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
DA_GO_MF_SemSim <- clusterSim(GO_genes[["LSS_DA"]], GO_genes[["ST_DA"]], ont= "MF",
                              organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
##### Test for similarity in combined GO terms sets (BP)
LSS_GO_BP_SemSim <- clusterSim(GO_genes[["LSS_DEGs"]], GO_genes[["LSS_DA"]], ont= "BP",
                               organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
ST_GO_BP_SemSim <- clusterSim(GO_genes[["ST_DEGs"]], GO_genes[["ST_DA"]], ont= "BP",
                              organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
DEGs_GO_BP_SemSim <- clusterSim(GO_genes[["LSS_DEGs"]], GO_genes[["ST_DEGs"]], ont= "BP",
                                organism= "human", measure= "Wang", drop= NULL, combine= "BMA")
DA_GO_BP_SemSim <- clusterSim(GO_genes[["LSS_DA"]], GO_genes[["ST_DA"]], ont= "BP",
                              organism= "human", measure= "Wang", drop= NULL, combine= "BMA")

##########################################################################################
##### ID GO & Pathway terms common to DEG & DA genesets (GO_MF)
LSS_vs_ST.DEGs.GO_MF.df <- read.delim("1701.LSS_vs_ST.DEGs.GO_MF", quote="\"'")
LSS_vs_ST.DEGs.GO_MF.df <- subset(LSS_vs_ST.DEGs.GO_MF.df[,c(1,2,3,8,9)], !duplicated(ID))
LSS.DEGs.GO_MF.df       <- subset(LSS_vs_ST.DEGs.GO_MF.df, Cluster == "LSS")
ST.DEGs.GO_MF.df        <- subset(LSS_vs_ST.DEGs.GO_MF.df, Cluster == "ST")

LSS_vs_ST.DA.GO_MF.df <- read.delim("1701.LSS_vs_ST.DA.GO_MF", quote="\"'")
LSS_vs_ST.DA.GO_MF.df <- subset(LSS_vs_ST.DA.GO_MF.df[,c(1,2,3,8,9)], !duplicated(ID))
LSS.DA.GO_MF.df       <- subset(LSS_vs_ST.DA.GO_MF.df, Cluster == "LSS")
ST.DA.GO_MF.df        <- subset(LSS_vs_ST.DA.GO_MF.df, Cluster == "ST")

## merge GO_MF terms to determine overlapping
LSS.DEG_and_DA_merged.GO_MF <- merge(LSS.DEGs.GO_MF.df, LSS.DA.GO_MF.df,
                                     by.x= "ID", by.y= "ID")
write.table(LSS.DEG_and_DA_merged.GO_MF, "LSS.DEG_and_DA_merged.GO_MF", sep="\t", row.names= F)
#
ST.DEG_and_DA_merged.GO_MF  <- merge(ST.DEGs.GO_MF.df,   ST.DA.GO_MF.df,
                                     by.x= "ID", by.y= "ID")
write.table(ST.DEG_and_DA_merged.GO_MF, "ST.DEG_and_DA_merged.GO_MF", sep="\t", row.names= F)
#
##### ID GO & Pathway terms common to DEG & DA genesets (GO_BP)
LSS_vs_ST.DEGs.GO_BP.df <- read.delim("1701.LSS_vs_ST.DEGs.GO_BP", quote="\"'")
LSS_vs_ST.DEGs.GO_BP.df <- subset(LSS_vs_ST.DEGs.GO_BP.df[,c(1,2,3,8,9)], !duplicated(ID))
LSS.DEGs.GO_BP.df       <- subset(LSS_vs_ST.DEGs.GO_BP.df, Cluster == "LSS")
ST.DEGs.GO_BP.df        <- subset(LSS_vs_ST.DEGs.GO_BP.df, Cluster == "ST")

LSS_vs_ST.DA.GO_BP.df <- read.delim("1701.LSS_vs_ST.DA.GO_BP", quote="\"'")
LSS_vs_ST.DA.GO_BP.df <- subset(LSS_vs_ST.DA.GO_BP.df[,c(1,2,3,8,9)], !duplicated(ID))
LSS.DA.GO_BP.df       <- subset(LSS_vs_ST.DA.GO_BP.df, Cluster == "LSS")
ST.DA.GO_BP.df        <- subset(LSS_vs_ST.DA.GO_BP.df, Cluster == "ST")

## merge GO_BP terms to determine overlapping
LSS.DEG_and_DA_merged.GO_BP <- merge(LSS.DEGs.GO_BP.df, LSS.DA.GO_BP.df,
                                     by.x= "ID", by.y= "ID")
write.table(LSS.DEG_and_DA_merged.GO_BP, "LSS.DEG_and_DA_merged.GO_BP", sep="\t", row.names= F)
#
ST.DEG_and_DA_merged.GO_BP  <- merge(ST.DEGs.GO_BP.df,   ST.DA.GO_BP.df,
                                     by.x= "ID", by.y= "ID")
write.table(ST.DEG_and_DA_merged.GO_BP, "ST.DEG_and_DA_merged.GO_BP", sep="\t", row.names= F)

