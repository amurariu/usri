library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#Human HIV (3) dataset
#####

raw_counts_hiv <- "https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_ASVs_table.tsv"
conds_hiv <-"https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_metadata.tsv"
hiv <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
hiv <- as.matrix(hiv)
hiv <- hiv+1
conditions_h <- read.table(file=conds_hiv, sep='\t', row.names = 1, header = T)
hiv.conds <- data.frame(conditions_h) 
transposed_hiv<- t(hiv)

# function
hiv.data.DESeq <- des.fun(data = transposed_hiv, nloop = 100, conditions = hiv.conds$comparison)
save(hiv.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/hiv.data.deseq.Rda") 