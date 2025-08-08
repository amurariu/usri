library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"

mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
mar <- mar[,-ncol(mar)]
conditions_m <- read.table(file=conds_mar, sep='\t', row.names = 1, header = T)
mar.conds <- data.frame(conditions_m) 

# filter to remove any ASV present in <10 % of samples
mar.filt <- mar[apply(mar, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

mar.data_4.aldex2 <- ald2.fun(data=as.matrix(mar.filt), conditions=mar.conds$comparison, nloop=100, gamma=0.4)
save(mar.data_4.aldex2, file="/Volumes/data2/andreea/ext_analysis/mar.data.aldex2_4.Rda")