devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')


#####
#c.diff(1) dataset
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"

cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)
cdiff.conds <- data.frame(conditions_c) 

cdiff.data_5.aldex3 <- ald3.fun(data=as.matrix(cdiff), conds=cdiff.conds$comparison, nloop=100, gamma=0.5)
save(cdiff.data_5.aldex3, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.aldex3_5.Rda")
