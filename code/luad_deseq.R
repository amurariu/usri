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

y_luad <- DGEList(counts=luad, group=factor(conditions_lu))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts 

# run deseq then save file
luad.data.DESeq <- des.fun(data = luad.data, nloop = 100, conditions = luad.conds$conditions_lu) 
save(luad.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/luad.data.deseq.Rda")