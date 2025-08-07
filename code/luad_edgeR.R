library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
#LUAD dataset
#####

raw_counts_luad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.normal-tumor.pair.rawCount.tsv"
conds_luad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.conditions.tsv"
luad <- read.table(file=raw_counts_luad, header=T, row.names=1, sep='\t')
luad.conds <- as.vector(unlist(read.table(file=conds_luad, sep='\t'))) 

#edgeR
y_luad <- DGEList(counts=luad, group=factor(luad.conds))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts 

# run edgeR then save file
luad.data.edgeR <- edg.fun(luad.data, luad.conds, 100)
save(luad.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/luad.data.edger.Rda")