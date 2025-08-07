library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
#c.diff(1) dataset
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"
cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)
cdiff.conds <- data.frame(conditions_c) 

# filter to remove any ASV present in <10 % of samples
cdiff.filt <- cdiff[apply(cdiff, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

# function
cdiff.data.edgeR <- edg.fun(data = as.matrix(cdiff.filt), conditions = cdiff.conds$comparison, nloop = 100)
save(cdiff.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.edger.Rda") 