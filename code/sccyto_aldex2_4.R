library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

source('code/ald2.fun.R')

#####
# sc_cyto dataset
#####

raw_counts_sccyto <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sc_cytoT_memB_counts.txt"
conds_sccyto <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sc_cytoT_memB_metadata.txt"

sccyto <- read.table(file=raw_counts_sccyto, header=T, row.names=1, sep='\t')
conditions_s <- read.table(file=conds_sccyto, sep='\t', row.names = 1, header = T)
sccyto.conds <- data.frame(conditions_s) 

# filter to remove any gene present in <10 % of samples
sccyto.filt <- sccyto[apply(sccyto, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

# run using gamma = 0.4 and save to .Rda
sccyto.data_4.aldex2 <- ald2.fun(data=as.matrix(sccyto.filt), conditions=sccyto.conds$comparison, nloop=100, gamma=0.4)
save(sccyto.data_4.aldex2, file="/Volumes/data2/andreea/ext_analysis/sccyto.data.aldex2_4.Rda")
