library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#marine sediment dataset
#####

raw_counts_mar <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_ASVs_table.tsv"
conds_mar <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sw_sed_detender_metadata.tsv"

mar <- read.table(file=raw_counts_mar, header=T, row.names=1, sep='\t')
conditions_m <- read.table(file=conds_mar, sep='\t', row.names = 1, header = T)

# remove taxonomy column
mar <- mar[,-ncol(mar)]

# filter to remove any ASV present in <10 % of samples
mar.filt <- mar[apply(mar, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

mar.data.limma <- lim.fun(data = as.matrix(mar.filt), conditions = conditions_m$comparison, nloop = 100, mean = 0, prop_null = 0.95)
save(mar.data.limma, file="/Volumes/data2/andreea/ext_analysis/mar.data.limma.Rda")
