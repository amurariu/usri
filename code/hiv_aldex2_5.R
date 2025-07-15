library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#Human HIV (3) dataset
#####

raw_counts_hiv <- "https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_ASVs_table.tsv"
conds_hiv <-"https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_metadata.tsv"

hiv <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
conditions_h <- read.table(file=conds_hiv, sep='\t', row.names = 1, header = T)
hiv.conds <- data.frame(conditions_h) 

hiv.data_5.aldex2 <- ald2.fun(data=as.matrix(hiv), conditions=hiv.conds$conditions_h, nloop=100, gamma=0.5)
save(hiv.data_5.aldex2, file="/Volumes/data2/andreea/ext_analysis/hiv.data.aldex2_5.Rda")