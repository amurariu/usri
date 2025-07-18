devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_ASVs_table.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_metadata.csv"

fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
conditions_f <- read.table(file=conds_fre, sep='\t', row.names = 1, header = T)

fre.data_1.aldex3 <- ald3.fun(data=as.matrix(fre), conds=conditions_f$comparison, nloop=100, gamma=0.1, prop_null = 0.95, mean = 0)
save(fre.data_1.aldex3, file="/Volumes/data2/andreea/ext_analysis/fre.data.aldex3_1.Rda")
