library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#PRAD dataset
#####

raw_counts_prad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.normal-tumor.pair.rawCount.tsv"
conds_prad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.conditions.tsv"
prad <- read.table(file=raw_counts_prad, header=T, row.names=1, sep='\t')
conditions_pr <- as.vector(unlist(read.table(file=conds_prad, sep='\t')))
prad.conds <- data.frame(conditions_pr) 

y_prad <- DGEList(counts=prad, group=factor(conditions_pr))
keep_prad <- filterByExpr(y_prad)
y_prad <- y_prad[keep_prad,keep.lib.sizes=FALSE]
prad.data <- y_prad$counts #filtered base dataset

#save file - gamma = 0.2
prad.data_2.aldex2 <- ald2.fun(data=prad.data, conditions=prad.conds$conditions_pr, nloop=100, gamma = 0.2)
save(prad.data_2.aldex2, file="/Volumes/data2/andreea/ext_analysis/prad.data.aldex2_2.Rda")
