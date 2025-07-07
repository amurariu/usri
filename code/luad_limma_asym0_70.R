library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#LUAD dataset
#####

raw_counts_luad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.normal-tumor.pair.rawCount.tsv"
conds_luad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.conditions.tsv"
luad <- read.table(file=raw_counts_luad, header=T, row.names=1, sep='\t')
luad.conds <- as.vector(unlist(read.table(file=conds_luad, sep='\t'))) 

y_luad <- DGEList(counts=luad, group=factor(luad.conds))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts #filtered base dataset

luad.data.limma.asym0.70 <- lim.fun(luad.data, luad.conds, 100, mean = 0, prop_null = 0.70)
save(luad.data.limma.asym0.70, file="/Volumes/data2/andreea/ext_analysis/luad.data.limma.asym0.prop70.Rda") 