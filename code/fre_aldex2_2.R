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

fre.data_2.aldex2 <- ald2.fun(data=as.matrix(fre), conditions=fre.conds$conditions_c, nloop=100, gamma=0.2)
save(fre.data_2.aldex2, file="../ext_analysis/fre.data.aldex2_2.Rda")