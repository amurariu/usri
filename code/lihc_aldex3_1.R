devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

#####
#LIHC dataset
#####

raw_counts_lihc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.normal-tumor.pair.rawCount.tsv"
conds_lihc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.conditions.tsv"

lihc <- read.table(file=raw_counts_lihc, header=T, row.names=1, sep='\t')
conditions_li <- as.vector(unlist(read.table(file=conds_lihc, sep='\t'))) 
lihc.conds <- data.frame(conditions_li) 


y_lihc <- DGEList(counts=lihc, group=factor(conditions_li))
keep_lihc <- filterByExpr(y_lihc)
y_lihc <- y_lihc[keep_lihc,keep.lib.sizes=FALSE]
lihc.data <- y_lihc$counts #filtered base dataset

# run at gamma = 0.1 and save
lihc.data_1.aldex3 <- ald3.fun(data=lihc.data, conds=lihc.conds$conditions_li, nloop=100, gamma=0.1)
save(lihc.data_1.aldex3, file="/Volumes/data2/andreea/ext_analysis/lihc.data.aldex3_1.Rda")