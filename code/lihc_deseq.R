library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#LIHC dataset
#####

raw_counts_lihc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.normal-tumor.pair.rawCount.tsv"
conds_lihc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.conditions.tsv"
lihc <- read.table(file=raw_counts_lihc, header=T, row.names=1, sep='\t')
conditions_li <- as.vector(unlist(read.table(file=conds_lihc, sep='\t'))) 
lihc.conds <- data.frame(conditions_li) 

y_lihc <- DGEList(counts=lihc, group=factor(conditions_li))
keep_lihc <- filterByExpr(y_lihc)
y_lihc <- y_lihc[keep_lihc,keep.lib.sizes=FALSE]
lihc.data <- y_lihc$counts 

#LIHC function
lihc.data.DESeq <- des.fun(data = lihc.data, nloop = 100, conditions = lihc.conds$conditions_li)
save(lihc.data.DESeq, file="../ext_analysis/lihc.data.deseq.Rda")
