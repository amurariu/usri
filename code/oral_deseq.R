library(seqgendiff, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
# oral dataset
#####

raw_counts_oral <- "https://raw.githubusercontent.com/amurariu/usri/main/data/oral_data.txt"
conds_oral <-"https://raw.githubusercontent.com/amurariu/usri/main/data/oral_meta.txt"

oral <- read.table(file=raw_counts_oral, header=T, row.names=1, sep='\t')
conditions_o <- read.table(file=conds_oral, sep='\t', row.names = 1, header = T)

# filter to remove any OTU present in <10 % of samples (dataset does not throw
# the same error regarding zeros in too many rows, so no need to add 1 to all
# counts)
oral.filt <- oral[apply(oral, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

# run DESeq2 and save
oral.data.DESeq <- des.fun(data = as.matrix(oral.filt), nloop = 100, conditions = conditions_o$group)
save(oral.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/oral.data.deseq.Rda") 
