library(seqgendiff, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#freshwater treatment dataset
#####

raw_counts_fre <- "https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_ASVs_table.tsv"
conds_fre <-"https://raw.githubusercontent.com/amurariu/usri/main/data/Ji_WTP_DS_metadata.csv"
fre <- read.table(file=raw_counts_fre, header=T, row.names=1, sep='\t')
fre <- as.matrix(fre)
conditions_f <- read.table(file=conds_fre, sep='\t', row.names = 1, header = T)

# filter to remove any ASV present in <10 % of samples, then add 1 to all read
# counts (common practice to avoid errors in DESeq2, but statistically incorrect)
fre.filt <- fre[apply(fre, 1, function(x){length(which(x != 0))/length(x)} >0.1),]
fre.filt <- fre.filt+1

# function
fre.data.DESeq <- des.fun(data = fre.filt, nloop = 100, conditions = conditions_f$comparison)
save(fre.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/fre.data.deseq.Rda") 