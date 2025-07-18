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
taxa <- mar$taxonomy
mar$taxonomy <- NULL
conditions_m <- read.table(file=conds_mar, sep='\t', header = T)

mar.data_5.aldex3 <- ald3.fun(data=as.matrix(mar), conds=conditions_m$comparison, nloop=100, gamma=0.5, prop_null = 0.95, mean = 0)
save(mar.data_5.aldex3, file="/Volumes/data2/andreea/ext_analysis/mar.data.aldex3_5.Rda")
