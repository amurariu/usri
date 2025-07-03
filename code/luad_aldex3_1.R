devtools::load_all('~/Documents/github/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#LUAD dataset
#####

raw_counts_luad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.normal-tumor.pair.rawCount.tsv"
conds_luad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.conditions.tsv"
luad <- read.table(file=raw_counts_luad, header=T, row.names=1, sep='\t')
conditions_lu <- as.vector(unlist(read.table(file=conds_luad, sep='\t'))) 
luad.conds <- data.frame(conditions_lu) 

y_luad <- DGEList(counts=luad, group=factor(conditions_lu))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts #filtered base dataset

luad.data_1.aldex3 <- ald3.fun(data=luad.data, conds=luad.conds$conditions_lu, nloop=100, gamma=0.1)
save(luad.data_1.aldex3, file="/Volumes/data2/andreea/ext_analysis/luad.data.aldex3_1.Rda")