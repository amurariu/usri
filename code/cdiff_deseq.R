library(seqgendiff, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/des.fun.R')

#####
#c.diff(1) dataset
#####

raw_counts_cdiff <- "https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_ASVs_table.tsv"
conds_cdiff <-"https://raw.githubusercontent.com/amurariu/usri/main/data/cdi_schubert_metadata.tsv"
cdiff <- read.table(file=raw_counts_cdiff, header=T, row.names=1, sep='\t')
cdiff <- as.matrix(cdiff)
conditions_c <- read.table(file=conds_cdiff, sep='\t', row.names = 1, header = T)

# filter to remove any ASV present in <10 % of samples, then add 1 to all read
# counts (common practice to avoid errors in DESeq2, but statistically incorrect)
cdiff.filt <- cdiff[apply(cdiff, 1, function(x){length(which(x != 0))/length(x)} >0.1),]
cdiff.filt <- cdiff.filt+1

# function
cdiff.data.DESeq <- des.fun(data = cdiff.filt, nloop = 100, conditions = conditions_c$comparison)
save(cdiff.data.DESeq, file="/Volumes/data2/andreea/ext_analysis/cdiff.data.deseq.Rda") 