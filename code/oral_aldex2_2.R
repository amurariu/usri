library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

source('code/ald2.fun.R')

#####
# oral dataset
#####

raw_counts_oral <- "https://raw.githubusercontent.com/amurariu/usri/main/data/oral_data.txt"
conds_oral <-"https://raw.githubusercontent.com/amurariu/usri/main/data/oral_meta.txt"

oral <- read.table(file=raw_counts_oral, header=T, row.names=1, sep='\t')
conditions_o <- read.table(file=conds_oral, sep='\t', row.names = 1, header = T)
oral.conds <- data.frame(conditions_o) 

# filter to remove any ASV present in <10 % of samples
oral.filt <- oral[apply(oral, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

oral.data_2.aldex2 <- ald2.fun(data=as.matrix(oral.filt), conditions=oral.conds$comparison, nloop=100, gamma=0.2)
save(oral.data_2.aldex2, file="/Volumes/data2/andreea/ext_analysis/oral.data.aldex2_2.Rda")