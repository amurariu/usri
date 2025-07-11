library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_ASVs_table.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_metadata.csv"
fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
conditions_f <- read.table(file=conds_fre, sep='\t', row.names = 1, header = T)
fre.conds <- data.frame(conditions_f) 

# function
fre.data.edgeR <- edg.fun(data = as.matrix(fre), conditions = fre.conds$comparison, nloop = 100)
save(fre.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/fre.data.edgeR.Rda") 