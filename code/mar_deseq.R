library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"
mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
taxa <- mar$taxonomy
mar$taxonomy <- NULL
conditions_m <- read.table(file=conds_mar, sep='\t', header = T)
transposed_mar<- t(mar)

# function
mar.data.DESeq <- des.fun(data = as.matrix(mar), nloop = 100, conditions = mar.conds$comparison)
save(mar.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/mar.data.deseq.Rda") 