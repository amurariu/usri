library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#BRCA dataset
#####

#pull the data and filter
raw_counts_brca<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts_brca, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=conds_brca, sep='\t')))
brca.conds <- data.frame(conditions_b) 

#edgeR
y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset

#save file - gamma = 0.1
brca.data_1.aldex2 <- ald2.fun(data=brca.data, conditions=brca.conds$conditions_b, nloop=100, gamma=0.1)
save(brca.data_1.aldex2, file="/Volumes/data2/andreea/ext_analysis/brca.data.aldex2_1.Rda")