setwd("/Users/andreeamurariu/Documents/github/usri/")

library(edgeR,quietly=T) 
library(DESeq2,quietly=T)
library(ALDEx2,quietly=T)
#add ALDEx3

readCount<-read.table(file="data/imm_GSE91061_raw_counts_GRCh38.p13_NCBI.tsv", header = T, skip=35, sep='\t', row.names = 1)
meta <- read.table(file="data/metadata.txt", header=F, row.names=1, sep='\t') #new file with accession numbers and patient IDs

conditions <- rep("Pre", 109)
grep("_On",meta)
conditions[grep("_On",meta)] <- "On"
conds <- data.frame(conditions)
conds$condition <- factor(conds$condition)