devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)

source('code/ald3.fun.R')

#####
# sc_cyto dataset
#####

raw_counts_sccyto <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sc_cytoT_memB_counts.txt"
conds_sccyto <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sc_cytoT_memB_metadata.txt"

sccyto <- read.table(file=raw_counts_sccyto, header=T, row.names=1, sep='\t')
conditions_s <- read.table(file=conds_sccyto, sep='\t', row.names = 1, header = T)
sccyto.conds <- data.frame(conditions_s) 

# filter to remove any ASV present in <10 % of samples
sccyto.filt <- sccyto[apply(sccyto, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

sccyto.data_1.aldex3 <- ald3.fun(data=as.matrix(sccyto.filt), conds=conditions_s$comparison, nloop=100, gamma=0.1, prop_null = 0.95, mean = 0)
save(sccyto.data_1.aldex3, file="/Volumes/data2/andreea/ext_analysis/sccyto.data.aldex3_1.Rda")
