library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#PRAD dataset
#####

raw_counts_prad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.normal-tumor.pair.rawCount.tsv"
conds_prad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.conditions.tsv"
prad <- read.table(file=raw_counts_prad, header=T, row.names=1, sep='\t')
prad.conds <- as.vector(unlist(read.table(file=conds_prad, sep='\t'))) 

y_prad <- DGEList(counts=prad, group=factor(prad.conds))
keep_prad <- filterByExpr(y_prad)
y_prad <- y_prad[keep_prad,keep.lib.sizes=FALSE]
prad.data <- y_prad$counts #filtered base dataset

prad.data.limma.asym1.80 <- lim.fun(prad.data, prad.conds, 100, mean = 1, prop_null = 0.80)
save(prad.data.limma.asym1.80, file="/Volumes/data2/andreea/ext_analysis/prad.data.limma.asym1.prop80.Rda") 
