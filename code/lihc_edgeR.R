library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
#LIHC dataset
#####

raw_counts_lihc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.normal-tumor.pair.rawCount.tsv"
conds_lihc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.conditions.tsv"
lihc <- read.table(file=raw_counts_lihc, header=T, row.names=1, sep='\t')
lihc.conds <- as.vector(unlist(read.table(file=conds_lihc, sep='\t'))) 

#edgeR
y_lihc <- DGEList(counts=lihc, group=factor(lihc.conds))
keep_lihc <- filterByExpr(y_lihc)
y_lihc <- y_lihc[keep_lihc,keep.lib.sizes=FALSE]
lihc.data <- y_lihc$counts 

# run edgeR and save
lihc.data.edgeR <- edg.fun(lihc.data, lihc.conds, 100)
save(lihc.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/lihc.data.edger.Rda")