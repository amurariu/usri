library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_ASVs_table.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_metadata.csv"

fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
conditions_f <- read.table(file=conds_fre, sep='\t', row.names = 1, header = T)
fre.conds <- data.frame(conditions_f) 

fre.data_0.aldex2 <- ald2.fun(data=as.matrix(fre), conditions=fre.conds$comparison, nloop=100, gamma=1e-3)
save(fre.data_0.aldex2, file="/Volumes/data2/andreea/ext_analysis/fre.data.aldex2_0.Rda")