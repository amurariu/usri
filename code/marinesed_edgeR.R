library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"

mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
conditions_m <- read.table(file=conds_mar, sep='\t', row.names = 1, header = T)
mar.conds <- data.frame(conditions_m) 

# function
mar.data.edgeR <- edg.fun(data = as.matrix(mar), conditions = mar.conds$comparison, nloop = 100)
save(mar.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/mar.data.edgeR.Rda") 