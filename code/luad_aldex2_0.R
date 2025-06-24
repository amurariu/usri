library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

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

#save file - gamma=1e-3
luad.data_0.aldex2 <- ald2.fun(data=luad.data, conditions=luad.conds$conditions_lu, nloop=100, gamma=1e-3)
save(luad.data_0.aldex2, file="../ext_analysis/luad.data.aldex2_0.Rda")