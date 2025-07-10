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
conditions_c <- as.vector(unlist(read.table(file=conds_cdiff, sep='\t')))
cdiff.conds <- data.frame(conditions_c) 


transposed_cdiff<- t(cdiff)

# function
cdiff.data.DESeq <- des.fun(data = cdiff.data, nloop = 100, conditions = cdiff.conds$conditions_c)
save(cdiff.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.deseq.Rda") 