# processing 16S rRNA amplicon data and running differential abundance scale 
# benchmarking - DESeq2

# Scott Dos Santos
# Last edited: 6th July 2026

# This code processes publicly available 16S rRNA amplicon sequencing data for
# thinning and randomisation. Data are taken from Bian et al. (2017). mSphere;
# 2 (5): e00327-17. This study looked at gut microbiome differences between a
# range of different age groups, from kindergarten-aged children to centenarians
# and compared OTUs present in the respective groups. This analysis will use a
# comparison of kindergarteners vs. middle-schoolers using DESeq2.

#################################### setup ####################################

library(DESeq2)
library(seqgendiff)

source("code/des.fun.R")
source("code/get_confusion.R")

# load and bind input data frames of 16S rRNA amplicon data
load("~/Documents/GitHub/tiyani-aldex3/data/kin.Rda") # private repo
load("~/Documents/GitHub/tiyani-aldex3/data/mid.Rda") # private repo

counts <- cbind(kin,mid)
ages <- data.frame(group = c(rep("K", 103), rep("M", 114)))

# filter OTUs with seqgendiff, choosing the 1000 most "median-expressed" genes
set.seed(123)
counts.sub <- select_counts(mat = as.matrix(counts), nsamp = ncol(counts), 
                            ngene = 1000,  gselect = "max")

# check for and remove zero sum rows
any(rowSums(counts.sub) == 0) # returns TRUE
counts.sub <- counts.sub[-which(rowSums(counts.sub) == 0),]
any(rowSums(counts.sub) == 0) # returns FALSE; removed 6 rows

############################ thinning & con matrix ############################

# run DESeq2 function for randomising and thinning
`16S.data.DESeq` <-  des.fun(data = counts.sub, conditions = ages$group,
                             nloop=100, prop_null = 0.05, mean = 0, std = 2)

obj <- "16S.data.DESeq"

# save output to .Rda
save(list = obj, file = "~/Documents/GitHub/usri/data/16S.data.deseq.Rda" )

# generate confusion matrix for the 100 iterations of thinned/randomised data
# for the corresponding gamma value and save
cm.16S.DESeq <- get_confusion(`16S.data.DESeq`, prog = "DESeq", FDR = 0.05)
save(cm.16S.DESeq, file = "~/Documents/GitHub/usri/analysis/confusionMats/cm.16S.deseq.Rda" )
