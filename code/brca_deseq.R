library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#BRCA dataset
#####

#pull the data and filter as above
raw_counts_brca<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts_brca, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=conds_brca, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- data.frame(conditions_b) 

y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset

# run deseq, then save
brca.data.DESeq <- des.fun(data = brca.data, nloop = 100, conditions = brca.conds$conditions_b) 
save(brca.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/brca.data.deseq.Rda")