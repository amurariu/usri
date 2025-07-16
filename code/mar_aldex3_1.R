devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"

mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
mar <- mar[,-ncol(mar)]
conditions_m <- as.vector(unlist(read.table(file=conds_mar, sep='\t')))
mar.conds <- data.frame(conditions_m) 

mar.data_1.aldex3 <- ald3.fun(data=mar, conds=mar.conds$comparison, nloop=100, gamma=0.1)
save(mar.data_1.aldex3, file="/Volumes/data2/andreea/ext_analysis/mar.data.aldex3_1.Rda")
