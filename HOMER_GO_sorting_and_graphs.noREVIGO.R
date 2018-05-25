########## Get gene ontologies, subset the important ones, make graphs
##### Read in files ####
##### Set title as the directory name
wd <- getwd()
Title <- sapply(strsplit(wd, split="\\/"), tail, 1L)
Title <- gsub('.HOMER_GO', '', Title)

## read all text files (various HOMER GO outputs) from CWD into a list
print(filenames <- list.files(pattern="*.txt"))
my_data <- lapply(filenames, function(files) {
  read.delim(files, header= T, sep = '\t')
})
## change the names of each df in list
names(my_data) <- gsub("\\.txt$", "", filenames)
## send all dataframes to GlobalEnv
list2env(my_data,.GlobalEnv)
print(testfiles <- names(my_data))



##### Create function to filter GO outputs & create graph ####
require(dplyr)
require(ggplot2)
require(export)
# df <- biological_process #test term
## Specify the number of terms to graph (default = 10)
FilterAndPlotGOterms <- function(df, dfName, termNum = 10) {
  TestTerm <- dfName
  #### Filter by Enrichment, Target.Genes.in.Term, Fraction.of.Targets.in.Term
  ## remove terms with Enrichment > 0.1
  df <- subset(df, Enrichment > 0.1)
  df <- subset(df, Target.Genes.in.Term > 1)
  ## remove terms who target >20% of tested terms
  number_tested_terms <- unique(df$Total.Target.Genes) * 0.2
  df <- subset(df, Target.Genes.in.Term < number_tested_terms)
  ## make p-val positive, sort by p-value
  df$logP <- df$logP * -1
  df <- arrange(df, desc(logP))
  ## subset df by desired number of terms to plot
  df <- df[1:termNum,] 
  ## set Terms as factors to prevent ABC order
  df$Term <- factor(df$Term, levels= df$Term)

  ## set plot features
  size <- element_text(size= 14) #font size on plot
  plotTitle <- paste(TestTerm, " (", unique(df$Total.Target.Genes), ")", sep='') 

  
  ## plot
  plota <- ggplot() + geom_bar(aes(y= logP, x= Term), data= df, stat="identity") +
    coord_flip() + ggtitle(plotTitle) + theme(plot.title= element_text(size= 14, face= "bold"),
                                          axis.text= size, legend.text= size,
                                          legend.title= size, axis.title= size) +
    geom_text(data=df, aes(x=Term, y=logP, label=as.factor(Target.Genes.in.Term)),hjust=-0.5)
  plot(plota)
  graph2ppt(file=paste(Title,".GO_and_Pathways.ppt",sep=''), width=10, height=5, append=T)
}

##### run function on all data.frames (GO Terms) ####
for (i in testfiles) {
  name      <- print(i)
  testframe <- get(i)
  FilterAndPlotGOterms(df = testframe, dfName = name, termNum = 10)
}
