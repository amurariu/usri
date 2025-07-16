library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_ASVs_table.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_metadata.csv"
fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
fre <- as.matrix(fre)
fre <- fre+1
conditions_f <- read.table(file=conds_fre, sep='\t', row.names = 1, header = T)
fre.conds <- data.frame(conditions_f) 
transposed_fre<- t(fre)

# function
fre.data.DESeq <- des.fun(data = transposed_fre, nloop = 100, conditions = fre.conds$comparison)
save(fre.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/fre.data.deseq.Rda") 