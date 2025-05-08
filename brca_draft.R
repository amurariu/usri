# you should not need this 
#setwd("/Users/andreeamurariu/Documents/github/usri/")
# load all your libraries first
library(edgeR,quietly=T) 
library(DESeq2,quietly=T)
library(seqgendiff, quietly=T) # fixed
library(ALDEx2, warn.conflicts=F)

##stopped editing here, start here

if(file.exists("analysis/thin_sim_brca.out.Rda")){
  load("analysis/thin_sim_brca.out.Rda") # file is data.out, list
} else {
  # fixed path, next two lines
  brca <- read.table('data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv', header=T, row.names=1, 
                     sep='\t')
  
  brca.conds <- as.vector(unlist(read.table('data/TCGA-BRCA.conditions.tsv', sep='\t')))
  conditions <- data.frame(brca.conds)
  
  
#   library(seqgendiff, warn.conflicts=F)
#   library(edgeR, warn.conflicts=F)
#   library(DESeq2, warn.conflicts=F)
  
  ###
  # edgeR functions
  
  y <- DGEList(counts=brca, group=factor(brca.conds))
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  
  # make the filtered base dataset
  brca.data <- y$counts
  
  brca.data.out <- list()