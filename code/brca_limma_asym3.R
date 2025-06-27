library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#BRCA dataset
#####

#pull the data and filter
raw_counts_brca<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts_brca, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=conds_brca, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- as.vector(conditions_b) #changed from conditions to brca.conds for consistency with PD1 dataset

y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset

brca.data.limma.asym3 <- lim.fun(brca.data, brca.conds, 100, mean = 3)
save(brca.data.limma.asym3, file="../ext_analysis/brca.data.limma.asym3.Rda") 