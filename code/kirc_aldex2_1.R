library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#KIRC
#####

raw_counts_kirc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.normal-tumor.pair.rawCount.tsv"
conds_kirc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.conditions.tsv"

kirc <- read.table(file=raw_counts_kirc, header=T, row.names=1, sep='\t')
conditions_k <- as.vector(unlist(read.table(file=conds_kirc, sep='\t'))) 
kirc.conds <- data.frame(conditions_k) #changed from conditions to brca.conds for consistency with PD1 dataset

y_kirc <- DGEList(counts=kirc, group=factor(conditions_k))
keep_kirc <- filterByExpr(y_kirc)
y_kirc <- y_kirc[keep_kirc,keep.lib.sizes=FALSE]
kirc.data <- y_kirc$counts #filtered base dataset

#save file - gamma = 0.1
kirc.data_1.aldex2 <- ald2.fun(data=kirc.data, conditions=kirc.conds$conditions_k, nloop=100, gamma=0.1)
save(kirc.data_1.aldex2, file="/Volumes/data2/andreea/ext_analysis/kirc.data.aldex2_1.Rda")