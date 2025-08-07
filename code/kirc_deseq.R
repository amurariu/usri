library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#KIRC
#####

raw_counts_kirc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.normal-tumor.pair.rawCount.tsv"
conds_kirc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.conditions.tsv"

kirc <- read.table(file=raw_counts_kirc, header=T, row.names=1, sep='\t')
conditions_k <- as.vector(unlist(read.table(file=conds_kirc, sep='\t'))) 
kirc.conds <- data.frame(conditions_k)

y_kirc <- DGEList(counts=kirc, group=factor(conditions_k))
keep_kirc <- filterByExpr(y_kirc)
y_kirc <- y_kirc[keep_kirc,keep.lib.sizes=FALSE]
kirc.data <- y_kirc$counts

# run deseq2 then save
kirc.data.DESeq <- des.fun(data = kirc.data, nloop = 100, conditions = kirc.conds$conditions_k)
save(kirc.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/kirc.data.deseq.Rda") 
