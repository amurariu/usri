library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#THCA dataset
#####

raw_counts_thca <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.normal-tumor.pair.rawCount.tsv"
conds_thca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.conditions.tsv"
thca <- read.table(file=raw_counts_thca, header=T, row.names=1, sep='\t')
thca.conds <- as.vector(unlist(read.table(file=conds_thca, sep='\t'))) 

y_thca <- DGEList(counts=thca, group=factor(thca.conds))
keep_thca <- filterByExpr(y_thca)
y_thca <- y_thca[keep_thca,keep.lib.sizes=FALSE]
thca.data <- y_thca$counts #filtered base dataset

thca.data.limma.asym1 <- lim.fun(thca.data, thca.conds, 100, mean = 1)
save(thca.data.limma.asym1, file="../ext_analysis/thca.data.limma.asym1.Rda") 