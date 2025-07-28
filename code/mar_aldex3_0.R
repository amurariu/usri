devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"

mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
taxa <- mar$taxonomy
mar$taxonomy <- NULL
conditions_m <- read.table(file=conds_mar, sep='\t', header = T)

# filter to remove any ASV present in <10 % of samples
mar.filt <- mar[apply(mar, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

mar.data_0.aldex3 <- ald3.fun(data=as.matrix(mar.filt), conds=conditions_m$comparison, nloop=100, gamma=1e-3, prop_null = 0.95, mean = 0)
save(mar.data_0.aldex3, file="/Volumes/data2/andreea/ext_analysis/mar.data.aldex3_0.Rda")
