library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#c.diff(1) datasetcd
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"

cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)
cdiff.conds <- data.frame(conditions_c) 

cdiff.data_2.aldex2 <- ald2.fun(data=as.matrix(cdiff), conditions=cdiff.conds$conditions_c, nloop=100, gamma=0.2)
save(cdiff.data_2.aldex2, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.aldex2_2.Rda")