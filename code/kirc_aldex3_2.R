devtools::load_all('~/Documents/github/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

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

kirc.data_2.aldex3 <- ald3.fun(data=kirc.data, conds=kirc.conds$conditions_k, nloop=100, gamma=0.2)
save(kirc.data_2.aldex3, file="/Volumes/data2/andreea/ext_analysis/kirc.data.aldex3_2.Rda")