library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

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
prad.data <- y_prad$counts 

#PRAD function
prad.data.DESeq <- des.fun(data = prad.data, nloop = 100, conditions = prad.conds$conditions_pr)
save(prad.data.DESeq, file="../ext_analysis/prad.data.deseq.Rda") 