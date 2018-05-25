##### Pathview ::  Visualize effect sizes on KEGG pathways ####

##### inputs ####
## need a matrix containing gene names (Symbol) & fold change information.
## need to know the KEGG pathway ID (http://www.genome.jp/kegg/pathway.html)

## load libraries #note that adding other libraries may interfere with
## the "select function in pathview command"
zzz <- c("pathview")
lapply(zzz, require, character.only=T)

##### read in gene expression data ####
geneTable <- read.delim("DMCs.PAH_vs_con.targetGenes.txt", quote="\"'")
## select relevant columns
geneTable <- geneTable[,c("Symbol", "mean_meth.diff")]
## switch the sign of methylation data (since it is repressive)
geneTable$methDiff <- geneTable$mean_meth.diff * -1
## make rows unique
geneTable          <- unique(geneTable)
## set Symbol as rownames
rownames(geneTable) <- geneTable$Symbol
## select relevant columns generate matrix (rownames = Symbols, col1 = effect size)
geneMatrix <- as.matrix(geneTable[,c("methDiff")])
rownames(geneMatrix) <- geneTable$Symbol

##### Set experiment title (output prefix) ####
Title <- "PAH_v_Con"

##### select pathway to test ####
pathwayID <- "hsa04020"
pathName  <- "Calcium_signaling"
## Wnt_signaling      ## "hsa04310" 
## Hippo_signaling    ## "hsa04390"
## Notch_signaling    ## "hsa04330"
## TGFb_signaling     ## "hsa04350"
## Calcium_signaling  ## "hsa04020"
## Pathway as PNG
pv.out <- pathview(gene.data= geneMatrix[,1], pathway.id= pathwayID, gene.idtype= "SYMBOL",
                   species= "hsa", out.suffix= paste(Title,pathName,sep="."), kegg.native= T,
                   low = (gene = "red"), high = (gene = "green"))



## Pathway as PDF
pv.out <- pathview(gene.dat = geneMatrix[,1], pathway.id= pathwayID, species= "hsa",
                   gene.idtype= "SYMBOL", out.suffix = paste(Title,pathName,sep=".")
                   ,kegg.native= F, sign.pos= "bottomleft")

