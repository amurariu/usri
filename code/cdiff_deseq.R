library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#c.diff(1) dataset
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"
cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
cdiff <- as.matrix(cdiff)
cdiff <- cdiff+1
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)
cdiff.conds <- data.frame(conditions_c) 
transposed_cdiff<- t(cdiff)

# function
cdiff.data.DESeq <- des.fun(data = transposed_cdiff, nloop = 100, conditions = cdiff.conds$comparison)
save(cdiff.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.deseq.Rda") 