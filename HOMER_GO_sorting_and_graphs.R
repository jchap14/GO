##### Source & load libraries ####
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("DESeq2","pasilla","DESeq","pathview"))
# install.packages(c("rJava","ReporteRs","ReporteRsjars","ggplot2","rtable","xtable","VennDiagram"))
# install.packages(c("taRifx","devtools","dplyr","data.table"))
# devtools::install_github('tomwenseleers/export',local=F)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R") 
zzz<-c("pheatmap","grid","gplots","ggplot2","export","devtools","DESeq2","pasilla","Biobase","EBSeq","dplyr","data.table",
       "genefilter","FactoMineR","VennDiagram","DOSE","ReactomePA","org.Hs.eg.db","clusterProfiler","pathview")
lapply(zzz, require, character.only=T)
data(geneList)

##### Set title ####
Title <- "DEGs.PAHvCon.FDR10"
enriched_terms <- paste(Title,'.enriched_terms.txt', sep = '')

##### Read in the functional_enrichment files from HOMER w/ header junk removed (made manually from HTML) ####
Terms.df <- read.delim(enriched_terms, quote="\"", header=T)

## all functional terms here
GO_Tree_terms <- as.vector(unique(Terms.df$GO.Tree))

## subset by 1 functional term at a time into a new df
GO_bp.df <- subset(Terms.df, Terms.df$GO.Tree == "biological process")
GO_bp.df <- GO_bp.df[,c(5,1,3,2,4,6:10)]
## write out df for use w Revigo online tool
write.table(GO_bp.df, paste(Title, ".RevigoIDs", sep=''), col.names= T, row.names=F, sep='\t')

#######################################################################################
### Run Revigo (online, 1st 2 columns w "Small Option"); then import the results, & sort/graph

##### Read in Revigo output & filter ####
ReviGO <- read.delim("REVIGO.csv", quote="\"", header=T, sep=',')
## remove NULL (redundant) terms
ReviGO <- subset(ReviGO, plot_X != "null")
## merge w genes in GO terms
ReviGO_w_terms <- merge(ReviGO, GO_bp.df, by.x= "term_ID" , by.y= "GO.ID")
## remove terms who target >20% of tested terms
number_tested_terms <- unique(ReviGO_w_terms$X..of.Target.Genes) * 0.2
number_tested_terms
## subset interesting columns
# ReviGO_w_terms <- ReviGO_w_terms[ReviGO_w_terms$X..of.Target.Genes.in.Term < number_tested_terms,]
ReviGO_w_terms <- ReviGO_w_terms[ReviGO_w_terms$X..of.Target.Genes.in.Term < 5,] #manual selection
ReviGO_w_terms2 <- ReviGO_w_terms[,c(1,2,7,15,16,17,20)]
## sort by p-value
ReviGO_w_terms2 <- arrange(ReviGO_w_terms2, log10.p.value)
write.table(ReviGO_w_terms2, paste(Title, ".RevigoReducedGOBP", sep=''), col.names= T, row.names=F, sep='\t')

##### Graph the results (plot more terms if going to remove manually) ####
termNum <- "25"
df <- ReviGO_w_terms2
size <- element_text(size= 14) #font size on plot

## make p-val positive, sort by p-value
df$log10.p.value <- df$log10.p.value * -1
df <- arrange(df, desc(log10.p.value))

##### Graph results (plot more terms if going to remove manually)
df$description <- factor(df$description, levels= df$description) #set X as factor preventing ABC order
plota <- ggplot() + geom_bar(aes(y= log10.p.value, x= description), data= df, stat="identity") +
  coord_flip() + ggtitle("GO BP") + theme(plot.title= element_text(size= 14, face= "bold"),
                                          axis.text= size, legend.text= size,
                                          legend.title= size, axis.title= size) +
  geom_text(data=df, aes(x=description, y=log10.p.value, label=as.factor(X..of.Target.Genes.in.Term)),hjust=-0.5)
plot(plota)
graph2ppt(file=paste(Title,".GO_and_Pathways.ppt",sep=''), width=10, height=9, append=T)

###############################################################################################
##### Subset & make graph for Reactome, KEGG, WikiPathways, BIOCYC, Pathway Interaction DB ####
Pathways.df <- subset(Terms.df, GO.Tree== "REACTOME pathways" | GO.Tree== "KEGG pathways" |
                        GO.Tree== "WikiPathways" | GO.Tree== "Pathway Interaction DB"| GO.Tree== "BIOCYC pathways") 
df <- Pathways.df[,c(5,1,3,2,4,6:10)]

## remove non-unique terms here by keeping version w/ most genes
df <- arrange(df, desc(X..of.Target.Genes.in.Term), Term)
df <- subset(df, !duplicated(Term))

## remove terms who target >20% of tested terms
number_tested_terms <- 10 #max(unique(df$X..of.Target.Genes)) * 0.2
df <- subset(df, X..of.Target.Genes.in.Term < number_tested_terms &
               X..of.Target.Genes.in.Term > 0) #specify min here
## calculate -log10pVal & subset/reorg columns
df$log10pVal <- log10(df$P.value) * -1
df <- df[,c(3,11,10,5,7,6,1,8,9)]
## sort by # of Target genes in Term
df <- arrange(df, desc(X..of.Target.Genes.in.Term))
write.table(df, paste(Title, ".pathways", sep=''), col.names= T, row.names=F, sep='\t')

##### Plot Pathways results ####
df <- df[1:31,]  #plot more terms if going to remove manually
df <- arrange(df, desc(log10pVal)) ## sort by p-value
df$Term <- factor(df$Term, levels= df$Term) #set X as factor preventing ABC order
plota <- ggplot() + geom_bar(aes(y= log10pVal, x= Term), data= df, stat="identity") +
  coord_flip() + ggtitle("Pathways") + theme(plot.title= element_text(size= 14, face= "bold"),
                                             axis.text= size, legend.text= size,
                                             legend.title= size, axis.title= size) +
  geom_text(data=df, aes(x=Term, y=log10pVal, label=as.factor(X..of.Target.Genes.in.Term)),hjust=-0.5)
plot(plota)
graph2ppt(file=paste(Title,".GO_and_Pathways.ppt",sep=''), width=10, height=9, append=T)

