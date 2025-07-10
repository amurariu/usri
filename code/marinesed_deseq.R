library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"
mar <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
conditions_m <- read.table(file=conds_mar, sep='\t', row.names = 1, header = T)
mar.conds <- data.frame(conditions_m) 
transposed_mar<- t(mar)


# function
mar.data.DESeq <- des.fun(data = transposed_mar, nloop = 100, conditions = mar.conds$conditions_m)
save(mar.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/mar.data.deseq.Rda") 