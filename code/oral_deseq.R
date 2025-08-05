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
oral.conds <- data.frame(conditions_o) 

# filter to remove any ASV present in <10 % of samples
oral.filt <- oral[apply(oral, 1, function(x){length(which(x != 0))/length(x)} >0.1),]
oral.filt <- oral.filt+1

# function
oral.data.DESeq <- des.fun(data = oral.filt, nloop = 100, conditions = conditions_o$comparison)
save(oral.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/oral.data.deseq.Rda") 