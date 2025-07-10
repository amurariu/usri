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
conditions_h <- as.vector(unlist(read.table(file=conds_hiv, sep='\t')))
hiv.conds <- data.frame(conditions_h) 



# function
hiv.data.DESeq <- des.fun(data = hiv.data, nloop = 100, conditions = hiv.conds$conditions_h)
save(hiv.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/hiv.data.deseq.Rda") 