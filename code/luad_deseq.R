library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#LUAD dataset
#####

raw_counts_luad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.normal-tumor.pair.rawCount.tsv"
conds_luad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.conditions.tsv"
luad <- read.table(file=raw_counts_luad, header=T, row.names=1, sep='\t')
conditions_lu <- as.vector(unlist(read.table(file=conds_luad, sep='\t'))) 
luad.conds <- data.frame(conditions_lu) 

y_luad <- DGEList(counts=luad, group=factor(luad.conds))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts 

#LUAD function
luad.data.DESeq <- des.fun(data = luad.data, nloop = 100, conditions = luad.conds$conditions_lu) 
save(luad.data.DESeq, file="../ext_analysis/luad.data.deseq.Rda")