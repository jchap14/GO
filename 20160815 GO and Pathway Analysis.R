##### Source & load libraries
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("DESeq2","pasilla","DESeq"))
# install.packages(c("rJava","ReporteRs","ReporteRsjars","ggplot2","rtable","xtable","VennDiagram"))
# install.packages(c("taRifx","devtools","dplyr","data.table"))
# devtools::install_github('tomwenseleers/export',local=F)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R") 
zzz<-c("pheatmap","grid","gplots","ggplot2","export","devtools","DESeq2","pasilla","Biobase","EBSeq","dplyr","data.table",
       "genefilter","FactoMineR","VennDiagram","DOSE","ReactomePA","org.Hs.eg.db","clusterProfiler","pathview")
lapply(zzz, require, character.only=T)

data(geneList)
##### Read in the DEGs
DEGs <- read.delim('iPSC-EC_vs_PAEC.overlapping_DEGs_in_all.avg.txt', quote="\"'")

##### Read in the name conversion file
name_convert <- read.delim('GRCh37_Ensblgene-Symbol-coords_NO-EXCEL.txt', quote="\"'")
rownames(name_convert) <- name_convert[,1]; name_convert[,1] = NULL

######import tables containing Ensembl geneIDs & adds EntrezIDs from a db
## PAEC
DEGs.PAEC = DEGs[DEGs$log2FC.iPSCEC_v_PAEC < 0, ]
DEGs.PAEC.Ensembl <- as.vector(DEGs.PAEC$Ensembl_nodot)
intersected.Ens_Entrez_PAEC <- clusterProfiler::bitr(DEGs.PAEC.Ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
PAEC <- as.vector(intersected.Ens_Entrez_PAEC[['ENTREZID']])
## iPSCEC
DEGs.iPSCEC = DEGs[DEGs$log2FC.iPSCEC_v_PAEC > 0, ]
DEGs.iPSCEC.Ensembl <- as.vector(DEGs.iPSCEC$Ensembl_nodot)
intersected.Ens_Entrez_iPSCEC <- clusterProfiler::bitr(DEGs.iPSCEC.Ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
iPSCEC <- as.vector(intersected.Ens_Entrez_iPSCEC[['ENTREZID']])
## COMBINED
COMBINED <- c(PAEC, iPSCEC)

DEG_cluster <- list(PAEC=c(PAEC), iPSCEC=c(iPSCEC)) ## for comparing seperately (Section 1)
COMBO_cluster <- list(COMBINED=c(COMBINED))          ## for combined analysis   (Section 2)

###### SECTION 1: Test *seperately* DEGs for GO enrichment using clusterProfiler::compareCluster #########
###### enrichGO Biological Process
cc_eGO_BP <- clusterProfiler::compareCluster(geneCluster= DEG_cluster, fun= "enrichGO",
                                             ont= "BP", OrgDb='org.Hs.eg.db',
                                             pAdjustMethod= "BH", readable=T,
                                             pvalueCutoff= 0.1, qvalueCutoff= 0.1)
#
plot(cc_eGO_BP, showCategory=10, title="eGO_BP", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#"dropGO" removes specified GO level
cc_eGO_BP_drop_level_1 <- clusterProfiler::dropGO(cc_eGO_BP, level= 1, term= NULL)
plot(cc_eGO_BP_drop_level_1, showCategory=10, title="eGO_BP_drop1", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_2 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_1, level= 2, term= NULL)
plot(cc_eGO_BP_drop_level_2, showCategory=10, title="eGO_BP_drop2", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_3 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_2, level= 3, term= NULL)
plot(cc_eGO_BP_drop_level_3, showCategory=10, title="eGO_BP_drop3", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_4 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_3, level= 4, term= NULL)
plot(cc_eGO_BP_drop_level_4, showCategory=10, title="eGO_BP_drop4", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_5 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_4, level= 5, term= NULL)
plot(cc_eGO_BP_drop_level_5, showCategory=10, title="eGO_BP_drop5", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
##Write out desired eGO_BP summaries
write.table(summary(cc_eGO_BP), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_BP", sep="\t", row.names=F)
write.table(summary(cc_eGO_BP_drop_level_1), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_BP_drop1.txt", sep="\t", row.names=F)
write.table(summary(cc_eGO_BP_drop_level_2), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_BP_drop2.txt", sep="\t", row.names=F)
write.table(summary(cc_eGO_BP_drop_level_3), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_BP_drop3.txt", sep="\t", row.names=F)
write.table(summary(cc_eGO_BP_drop_level_4), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_BP_drop4.txt", sep="\t", row.names=F)
write.table(summary(cc_eGO_BP_drop_level_5), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_BP_drop5.txt", sep="\t", row.names=F)

##### enrichGO Molecular Function
cc_eGO_MF <- clusterProfiler::compareCluster(geneCluster= DEG_cluster, fun= "enrichGO",
                                             ont= "MF", OrgDb='org.Hs.eg.db',
                                             pAdjustMethod= "BH", readable=T,
                                             pvalueCutoff= 0.1, qvalueCutoff= 0.1)
## GENERATE PLOTS
plot(cc_eGO_MF, showCategory=10, title="eGO_MF", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_MF.ppt", width=9, height=7, append=T)
dev.off()
#"dropGO" removes specified GO level
cc_eGO_MF_drop_level_1 <- clusterProfiler::dropGO(cc_eGO_MF, level= 1, term= NULL)
plot(cc_eGO_MF_drop_level_1, showCategory=10, title="eGO_MF_drop1", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_MF.ppt", width=9, height=7, append=T)
dev.off()
#
cc_eGO_MF_drop_level_2 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_1, level= 2, term= NULL)
plot(cc_eGO_MF_drop_level_2, showCategory=10, title="eGO_MF_drop2", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_MF.ppt", width=9, height=7, append=T)
dev.off()
#
cc_eGO_MF_drop_level_3 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_2, level= 3, term= NULL)
plot(cc_eGO_MF_drop_level_3, showCategory=10, title="eGO_MF_drop3", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_MF.ppt", width=9, height=7, append=T)
dev.off()
#
cc_eGO_MF_drop_level_4 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_3, level= 4, term= NULL)
plot(cc_eGO_MF_drop_level_4, showCategory=10, title="eGO_MF_drop4", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_MF.ppt", width=9, height=7, append=T)
dev.off()
#
cc_eGO_MF_drop_level_5 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_4, level= 5, term= NULL)
plot(cc_eGO_MF_drop_level_5, showCategory=10, title="eGO_MF_drop5", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.GO_MF.ppt", width=9, height=7, append=T)
dev.off()

##Write out desired eGO_MF summaries
# write.table(summary(cc_eGO_MF), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_MF", sep="\t", row.names=F)
# write.table(summary(cc_eGO_MF_drop_level_1), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_MF_drop1.txt", sep="\t", row.names=F)
# write.table(summary(cc_eGO_MF_drop_level_2), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_MF_drop2.txt", sep="\t", row.names=F)
# write.table(summary(cc_eGO_MF_drop_level_3), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_MF_drop3.txt", sep="\t", row.names=F)
write.table(summary(cc_eGO_MF_drop_level_4), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_MF_drop4.txt", sep="\t", row.names=F)
# write.table(summary(cc_eGO_MF_drop_level_5), "iPSC-EC_vs_PAEC.overlapping_DEGs.eGO_MF_drop5.txt", sep="\t", row.names=F)

##### enrichDO
cc_eDO <- clusterProfiler::compareCluster(geneCluster= DEG_cluster, fun= "enrichDO",
                                          pvalueCutoff= 0.1, readable=T)
plot(cc_eDO, showCategory=10, title="Disease Ontology", by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.pathways.ppt", width=10, height=6, append=T)
dev.off()
write.table(summary(cc_eDO), "iPSC-EC_vs_PAEC.overlapping_DEGs.DO.txt", sep="\t", row.names=F)

#####  Reactome
cc_reactome <- clusterProfiler::compareCluster(geneCluster= DEG_cluster, fun= "ReactomePA::enrichPathway",
                                               organism= "human", pvalueCutoff= 0.1, readable=T)
plot(cc_reactome, showCategory=8, by= "count", font.size= 14, title="Reactome")
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.pathways.ppt", width=6, height=3, append=T)
dev.off()
write.table(summary(cc_reactome), "iPSC-EC_vs_PAEC.overlapping_DEGs.reactome.txt", sep="\t", row.names=F)

#####  KEGG
cc_KEGG <- clusterProfiler::compareCluster(geneCluster= DEG_cluster, fun= "enrichKEGG", 
                                           use_internal_data= F, organism= "human",
                                           pvalueCutoff= 0.1)
cc_KEGG_summary <- (summary(cc_KEGG))
write.table(cc_KEGG_summary, "iPSC-EC_vs_PAEC.overlapping_DEGs.KEGG.txt", sep="\t", row.names=F)
plot(cc_KEGG, showCategory=10, by= "count", font.size= 14)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.pathways.ppt", width=10, height=6, append=T)
dev.off()

###### SECTION 2: Test *combined* DEGs for GO enrichment using clusterProfiler::compareCluster #########
###### enrichGO Biological Process
cc_eGO_BP <- clusterProfiler::compareCluster(geneCluster = COMBO_cluster, fun = "enrichGO",
                                             ont= "BP", OrgDb='org.Hs.eg.db',
                                             pAdjustMethod= "BH",
                                             pvalueCutoff= 0.1, qvalueCutoff= 0.1)
cc_eGO_BP_summary <- (summary(cc_eGO_BP))
write.table(cc_eGO_BP_summary, "iPSC-EC_vs_PAEC.combined_overlapping_DEGs.eGO_BP", sep="\t", row.names=F)
#
plot(cc_eGO_BP, showCategory=10, title="eGO_BP")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#"dropGO" removes specified GO level
cc_eGO_BP_drop_level_1 <- clusterProfiler::dropGO(cc_eGO_BP, level = 1, term= NULL)
plot(cc_eGO_BP_drop_level_1, showCategory=10, title="eGO_BP_drop1")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_2 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_1, level = 2, term= NULL)
plot(cc_eGO_BP_drop_level_2, showCategory=10, title="eGO_BP_drop2")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_3 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_2, level = 3, term= NULL)
plot(cc_eGO_BP_drop_level_3, showCategory=10, title="eGO_BP_drop3")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_4 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_3, level = 4, term= NULL)
plot(cc_eGO_BP_drop_level_4, showCategory=10, title="eGO_BP_drop4")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_BP_drop_level_5 <- clusterProfiler::dropGO(cc_eGO_BP_drop_level_4, level = 5, term= NULL)
plot(cc_eGO_BP_drop_level_5, showCategory=10, title="eGO_BP_drop5")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_BP.ppt", width=10, height=6, append=T)
dev.off()

##### enrichGO Molecular Function
cc_eGO_MF <- clusterProfiler::compareCluster(geneCluster = COMBO_cluster, fun = "enrichGO",
                                             ont= "MF", OrgDb='org.Hs.eg.db',
                                             pAdjustMethod = "BH",
                                             pvalueCutoff = 0.1, qvalueCutoff= 0.1)
## write out summary
cc_eGO_MF_summary <- (summary(cc_eGO_MF))
write.table(cc_eGO_MF_summary, "iPSC-EC_vs_PAEC.combined_overlapping_DEGs.eGO_MF", sep="\t", row.names=F)
## export dotplots to PPT
plot(cc_eGO_MF, showCategory=10, title="eGO_MF")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_MF.ppt", width=10, height=6, append=T)
dev.off()
## "dropGO" removes specified GO level
cc_eGO_MF_drop_level_1 <- clusterProfiler::dropGO(cc_eGO_MF, level = 1, term= NULL)
plot(cc_eGO_MF_drop_level_1, showCategory=10, title="eGO_MF_drop1")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_MF.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_MF_drop_level_2 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_1, level = 2, term= NULL)
plot(cc_eGO_MF_drop_level_2, showCategory=10, title="eGO_MF_drop2")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_MF.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_MF_drop_level_3 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_2, level = 3, term= NULL)
plot(cc_eGO_MF_drop_level_3, showCategory=10, title="eGO_MF_drop3")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_MF.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_MF_drop_level_4 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_3, level = 4, term= NULL)
plot(cc_eGO_MF_drop_level_4, showCategory=10, title="eGO_MF_drop4")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_MF.ppt", width=10, height=6, append=T)
dev.off()
#
cc_eGO_MF_drop_level_5 <- clusterProfiler::dropGO(cc_eGO_MF_drop_level_4, level = 5, term= NULL)
plot(cc_eGO_MF_drop_level_5, showCategory=10, title="eGO_MF_drop5")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.GO_MF.ppt", width=10, height=6, append=T)
dev.off()
#

##### enrichDO
cc_eDO <- clusterProfiler::compareCluster(geneCluster= COMBO_cluster, fun= "enrichDO",
                                          pvalueCutoff= 0.1)
cc_eDO_summary <- (summary(cc_eDO))
write.table(cc_eDO_summary, "iPSC-EC_vs_PAEC.overlapping_DEGs.DO.txt", sep="\t", row.names=F)
plot(cc_eDO, showCategory=10, title="Disease Ontology")
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.pathways.ppt", width=10, height=6, append=T)
dev.off()

#####  Reactome
cc_reactome <- clusterProfiler::compareCluster(geneCluster= COMBO_cluster, fun= "ReactomePA::enrichPathway",
                                               organism= "human", pvalueCutoff= 0.1)
cc_reactome_summary <- (summary(cc_reactome))
write.table(cc_reactome_summary, "iPSC-EC_vs_PAEC.overlapping_DEGs.reactome.txt", sep="\t", row.names=F)
plot(cc_reactome, showCategory=10)
graph2ppt(file="iPSC-EC_vs_PAEC.combined_overlapping_DEGs.pathways.ppt", width=12, height=4, append=T)
dev.off()

#####  KEGG
cc_KEGG <- clusterProfiler::compareCluster(geneCluster= COMBO_cluster, fun= "enrichKEGG", 
                                           use_internal_data= F, organism = "human",
                                           pvalueCutoff= 0.1)
cc_KEGG_summary <- (summary(cc_KEGG))
write.table(cc_KEGG_summary, "iPSC-EC_vs_PAEC.combined_overlapping_DEGs.KEGG.txt", sep="\t", row.names=F)
plot(cc_KEGG, showCategory=10)
graph2ppt(file="iPSC-EC_vs_PAEC.overlapping_DEGs.pathways.ppt", width=10, height=6, append=T)
dev.off()
