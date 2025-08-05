devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)

source('code/ald3.fun.R')

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

oral.data_0.aldex3 <- ald3.fun(data=as.matrix(oral.filt), conds=conditions_s$comparison, nloop=100, gamma=1e-3, prop_null = 0.95, mean = 0)
save(oral.data_0.aldex3, file="/Volumes/data2/andreea/ext_analysis/oral.data.aldex3_0.Rda")
