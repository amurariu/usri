devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#Human HIV (3) dataset
#####

raw_counts_hiv <- "https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_ASVs_table.tsv"
conds_hiv <-"https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_metadata.tsv"

hiv <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
conditions_h <- read.table(file=conds_hiv, sep='\t', row.names = 1, header = T)
hiv.conds <- data.frame(conditions_h) 

hiv.data_1.aldex3 <- ald3.fun(data=as.matrix(hiv), conds=hiv.conds$comparison, nloop=100, gamma=0.1)
save(hiv.data_1.aldex3, file="/Volumes/data2/andreea/ext_analysis/hiv.data.aldex3_1.Rda")
