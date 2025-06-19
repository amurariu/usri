library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
# PD1 immunotherapy dataset
####

# pull the raw count data from github
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'

immuno <- read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')

# make conditions vector for PD1 dataset using metadata info (Pre vs. On)
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p)

# edgeR conditions for initial filtering
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset

#Immuno function
immuno.data.DESeq <- des.fun(data = immuno.data, nloop = 100, conditions = immuno.conds$conditions_p)
save(immuno.data.DESeq, file="../ext_analysis/immuno.data.deseq.Rda") 

