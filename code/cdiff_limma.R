library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#c.diff(1) dataset
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"
cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)
cdiff.conds <- as.vector(unlist(conditions_c))
transposed_cdiff<- t(cdiff)

cdiff.data.limma <- lim.fun(data = transposed_cdiff, conditions = cdiff.conds, nloop = 100, mean = 0, prop_null = 0.95)
save(cdiff.data.limma, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.limma.Rda") 

