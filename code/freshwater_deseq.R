library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/name.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/name.tsv"
fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
conditions_f <- as.vector(unlist(read.table(file=conds_fre, sep='\t')))
fre.conds <- data.frame(conditions_f) 

y_fre <- DGEList(counts=fre, group=factor(conditions_h))
keep_fre <- filterByExpr(y_fre)
y_fre <- y_fre[keep_fre,keep.lib.sizes=FALSE]
fre.data <- y_fre$counts 

# function
fre.data.DESeq <- des.fun(data = fre.data, nloop = 100, conditions = fre.conds$conditions_f)
save(fre.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/fre.data.deseq.Rda") 