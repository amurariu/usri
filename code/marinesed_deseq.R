library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/name.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/name.tsv"
mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
conditions_m <- as.vector(unlist(read.table(file=conds_mar, sep='\t')))
mar.conds <- data.frame(conditions_m) 

y_mar <- DGEList(counts=mar, group=factor(conditions_m))
keep_mar <- filterByExpr(y_mar)
y_mar <- y_mar[keep_mar,keep.lib.sizes=FALSE]
mar.data <- y_mar$counts 

# function
mar.data.DESeq <- des.fun(data = mar.data, nloop = 100, conditions = mar.conds$conditions_m)
save(mar.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/mar.data.deseq.Rda") 