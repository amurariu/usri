# processing 16S rRNA amplicon data and running differential abundance scale 
# benchmarking - edgeR

# Scott Dos Santos
# Last edited: 15th July 2026

# This code processes publicly available 16S rRNA amplicon sequencing data for
# thinning and randomisation. Data are taken from Bian et al. (2017). mSphere;
# 2 (5): e00327-17. This study looked at gut microbiome differences between a
# range of different age groups, from kindergarten-aged children to centenarians
# and compared OTUs present in the respective groups. This analysis will use a
# comparison of kindergarteners vs. middle-schoolers using edgeR.

#################################### setup ####################################

library(edgeR)
library(seqgendiff)

source("code/edg.fun.R")
source("code/get_confusion.R")

# load input data frames of 16S rRNA amplicon data
rda <- paste0("https://github.com/amurariu/usri/raw/refs/heads/main/data/tianyi_",
              c("mid", "kin"), "_ft.Rda")

for(i in rda){
  tf <- tempfile(fileext = ".Rda")
  download.file(url = i, destfile = tf, mode = "wb")
  load(tf)
  unlink(tf)
  rm(tf)
}

# bind input data and make metadata objects
counts <- cbind(kin,mid)
ages <- data.frame(group = c(rep("K", 103), rep("M", 114)))

# filter OTUs with seqgendiff, choosing the 1000 most "median-expressed" genes
set.seed(123)
counts.sub <- select_counts(mat = as.matrix(counts), nsamp = ncol(counts), 
                            ngene = 400,  gselect = "max")

# check for and remove zero sum rows
any(rowSums(counts.sub) == 0) # returns FALSE
# counts.sub <- counts.sub[-which(rowSums(counts.sub) == 0),]
# any(rowSums(counts.sub) == 0) # returns FALSE

############################ thinning & con matrix ############################

# run edgeR function for randomising and thinning
`16S.data.edgeR` <-  edg.fun(data = counts.sub, conditions = ages$group,
                             nloop=100, prop_null = 0.95, mean = 0, std = 4)

obj <- "16S.data.edgeR"

# save output to .Rda
save(list = obj, file = "~/Documents/GitHub/ext_analysis/16S.data.edger.Rda" )

# generate confusion matrix for the 100 iterations of thinned/randomised data
# for the corresponding gamma value and save
cm.16S.edgeR <- get_confusion(`16S.data.edgeR`, prog = "edgeR", FDR = 0.05)
save(cm.16S.edgeR, file = "~/Documents/GitHub/usri/analysis/confusionMats/cm.16S.edger.Rda" )
