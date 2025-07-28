library(seqgendiff, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#Human HIV (3) dataset
#####

raw_counts_hiv <- "https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_ASVs_table.tsv"
conds_hiv <-"https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_metadata.tsv"
hiv <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
hiv <- as.matrix(hiv)
conditions_h <- read.table(file=conds_hiv, sep='\t', row.names = 1, header = T)

# filter to remove any ASV present in <10 % of samples, then add 1 to all read
# counts (common practice to avoid errors in DESeq2, but statistically incorrect)
hiv.filt <- hiv[apply(hiv, 1, function(x){length(which(x != 0))/length(x)} >0.1),]
hiv.filt <- hiv.filt +1

# function
hiv.data.DESeq <- des.fun(data = hiv.filt, nloop = 100, conditions = conditions_h$comparison)
save(hiv.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/hiv.data.deseq.Rda") 