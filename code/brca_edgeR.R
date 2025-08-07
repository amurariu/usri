library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/edg.fun.R')

#####
#BRCA dataset
#####

#pull the data and filter
raw_counts_brca <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"
brca <- read.table(file=raw_counts_brca, header=T, row.names=1, sep='\t')
brca.conds <- as.vector(unlist(read.table(file=conds_brca, sep='\t'))) 

#edgeR
y_brca <- DGEList(counts=brca, group=factor(brca.conds))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts 

# run edgeR, then save
brca.data.edgeR <- edg.fun(brca.data, brca.conds, 100)
save(brca.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/brca.data.edger.Rda")