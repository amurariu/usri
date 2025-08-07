library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/edg.fun.R')

#####
# sc_cyto dataset
#####

raw_counts_sccyto <- "https://raw.githubusercontent.com/amurariu/usri/main/data/sc_cytoT_memB_counts.txt"
conds_sccyto <-"https://raw.githubusercontent.com/amurariu/usri/main/data/sc_cytoT_memB_metadata.txt"

sccyto <- read.table(file=raw_counts_sccyto, header=T, row.names=1, sep='\t')
conditions_s <- read.table(file=conds_sccyto, sep='\t', row.names = 1, header = T)

# filter to remove any gene present in <10 % of samples
sccyto.filt <- sccyto[apply(sccyto, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

# convert to matrix
sccyto.filt <- as.matrix(sccyto.filt)

# run edgeR then save
sccyto.data.edgeR <- edg.fun(data = sccyto.filt, conditions = conditions_s$group, nloop = 100)
save(sccyto.data.edgeR, file="/Volumes/data2/andreea/ext_analysis/sccyto.data.edger.Rda") 