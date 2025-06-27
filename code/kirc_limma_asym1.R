library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#KIRC
#####

raw_counts_kirc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.normal-tumor.pair.rawCount.tsv"
conds_kirc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.conditions.tsv"

kirc <- read.table(file=raw_counts_kirc, header=T, row.names=1, sep='\t')
kirc.conds <- as.vector(unlist(read.table(file=conds_kirc, sep='\t'))) 

y_kirc <- DGEList(counts=kirc, group=factor(kirc.conds))
keep_kirc <- filterByExpr(y_kirc)
y_kirc <- y_kirc[keep_kirc,keep.lib.sizes=FALSE]
kirc.data <- y_kirc$counts #filtered base dataset

kirc.data.limma.asym1 <- lim.fun(kirc.data, kirc.conds, 100, mean = 1)
save(kirc.data.limma.asym1, file="../ext_analysis/kirc.data.limma.asym1.Rda") 
