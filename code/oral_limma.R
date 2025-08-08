library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
# oral dataset
#####

raw_counts_oral <- "https://raw.githubusercontent.com/amurariu/usri/main/data/oral_data.txt"
conds_oral <-"https://raw.githubusercontent.com/amurariu/usri/main/data/oral_meta.txt"

oral <- read.table(file=raw_counts_oral, header=T, row.names=1, sep='\t')
conditions_o <- read.table(file=conds_oral, sep='\t', row.names = 1, header = T)

# filter to remove any OTU present in <10 % of samples
oral.filt <- oral[apply(oral, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

# run limma, then save
oral.data.limma <- lim.fun(data = as.matrix(oral.filt), conditions = conditions_o$group, nloop = 100, mean = 0, prop_null = 0.95)
save(oral.data.limma, file="/Volumes/data2/andreea/ext_analysis/oral.data.limma.Rda")
