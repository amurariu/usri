devtools::load_all('~/Documents/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

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

#save file - gamma = 1e-3
thca.data_0.aldex3 <- ald3.fun(data=thca.data, conds=thca.conds$conditions_t, nloop=100, gamma=1e-3)
save(thca.data_0.aldex3, file="/Volumes/data2/andreea/ext_analysis/thca.data.aldex3_0.Rda")
