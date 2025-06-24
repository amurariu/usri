library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/edg.fun.R')

#####
##PD1 Immuno dataset
#####

#loading and filtering of data
raw_counts_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno <-read.table(file=raw_counts_immuno, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta_immuno, header=F, row.names=1, sep='\t')
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- as.vector(conditions_p)

#edgeR 
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset

immuno.data.edgeR <-edg.fun(immuno.data, immuno.conds,100)
save(immuno.data.edgeR, file="../ext_analysis/immuno.data.edger.Rda")