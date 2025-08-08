devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)

source('code/ald3.fun.R')


#####
#c.diff(1) dataset
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"

cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)

# filter to remove any ASV present in <10 % of samples
cdiff.filt <- cdiff[apply(cdiff, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

cdiff.data_3.aldex3 <- ald3.fun(data=as.matrix(cdiff.filt), conds=conditions_c$comparison, nloop=100, gamma=0.3, prop_null = 0.95, mean = 0)
save(cdiff.data_3.aldex3, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.aldex3_3.Rda")
