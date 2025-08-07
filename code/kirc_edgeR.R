library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/edg.fun.R')

#####
#KIRC dataset
#####

raw_counts_kirc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.normal-tumor.pair.rawCount.tsv"
conds_kirc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.conditions.tsv"
kirc <- read.table(file=raw_counts_kirc, header=T, row.names=1, sep='\t')
kirc.conds <- as.vector(unlist(read.table(file=conds_kirc, sep='\t'))) 

#edgeR
y_kirc <- DGEList(counts=kirc, group=factor(kirc.conds))
keep_kirc <- filterByExpr(y_kirc)
y_kirc <- y_kirc[keep_kirc,keep.lib.sizes=FALSE]
kirc.data <- y_kirc$counts 

# run edgeR then save
kirc.data.edgeR <- edg.fun(kirc.data, kirc.conds, 100)
save(kirc.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/kirc.data.edger.Rda")