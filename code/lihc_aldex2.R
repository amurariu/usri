library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

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

#save file - gamma=1e-3
lihc.data_0.aldex2 <- ald2.fun(data=lihc.data, conditions=lihc.conds, nloop=100, gamma=1e-3)
save(lihc.data_0.aldex2, file="../ext_analysis/lihc.data.aldex2_0.out.Rda")

#save file - gamma=0.2
lihc.data_2.aldex2 <- ald2.fun(data=lihc.data, conditions=lihc.conds, nloop=100, gamma = 0.2)
save(lihc.data_2.aldex2, file="../ext_analysis/lihc.data.aldex2_2.out.Rda")

#save file - gamma=0.5
lihc.data_5.aldex2 <- ald2.fun(data=lihc.data, conditions=lihc.conds, nloop=100, gamma = 0.5)
save(lihc.data_5.aldex2, file="../ext_analysis/lihc.data.aldex2_5.out.Rda")