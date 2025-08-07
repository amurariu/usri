library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
#THCA dataset
#####

raw_counts_thca <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.normal-tumor.pair.rawCount.tsv"
conds_thca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.conditions.tsv"
thca <- read.table(file=raw_counts_thca, header=T, row.names=1, sep='\t')
thca.conds <- as.vector(unlist(read.table(file=conds_thca, sep='\t'))) 

#edgeR
y_thca <- DGEList(counts=thca, group=factor(thca.conds))
keep_thca <- filterByExpr(y_thca)
y_thca <- y_thca[keep_thca,keep.lib.sizes=FALSE]
thca.data <- y_thca$counts 

# run edgeR then save
thca.data.edgeR <- edg.fun(thca.data, thca.conds, 100)
save(thca.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/thca.data.edger.Rda")