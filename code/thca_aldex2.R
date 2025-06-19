library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#THCA dataset
#####

raw_counts_thca <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.normal-tumor.pair.rawCount.tsv"
conds_thca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.conditions.tsv"
thca <- read.table(file=raw_counts_thca, header=T, row.names=1, sep='\t')
conditions_t <- as.vector(unlist(read.table(file=conds_thca, sep='\t')))
thca.conds <- data.frame(conditions_t) 

y_thca <- DGEList(counts=thca, group=factor(conditions_t))
keep_thca <- filterByExpr(y_thca)
y_thca <- y_thca[keep_thca,keep.lib.sizes=FALSE]
thca.data <- y_thca$counts #filtered base dataset

#save file - gamma=1e-3
thca.data_0.aldex2 <- ald2.fun(data=thca.data, conditions=thca.conds, nloop=100, gamma=1e-3)
save(thca.data_0.aldex2, file="../ext_analysis/thca.data.aldex2_0.out.Rda")

#save file - gamma=0.2
thca.data_2.aldex2 <- ald2.fun(data=thca.data, conditions=thca.conds, nloop=100, gamma = 0.2)
save(thca.data_2.aldex2, file="../ext_analysis/thca.data.aldex2_2.out.Rda")

#save file - gamma=0.5
thca.data_5.aldex2 <- ald2.fun(data=thca.data, conditions=thca.conds, nloop=100, gamma = 0.5)
save(thca.data_5.aldex2, file="../ext_analysis/thca.data.aldex2_5.out.Rda")
