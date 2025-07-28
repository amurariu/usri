devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_ASVs_table.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_metadata.csv"

fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
conditions_f <- read.table(file=conds_fre, sep='\t', row.names = 1, header = T)

# filter to remove any ASV present in <10 % of samples
fre.filt <- fre[apply(fre, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

fre.data_0.aldex3 <- ald3.fun(data=as.matrix(fre.filt), conds=conditions_f$comparison, nloop=100, gamma=1e-3, prop_null = 0.95, mean = 0)
save(fre.data_0.aldex3, file="/Volumes/data2/andreea/ext_analysis/fre.data.aldex3_0.Rda")
