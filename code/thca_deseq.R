library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

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
thca.data <- y_thca$counts 

#THCA function
thca.data.DESeq <- des.fun(data = thca.data, nloop = 100, conditions = thca.conds$conditions_t)
save(thca.data.DESeq, file="../ext_analysis/thca.data.deseq.Rda")