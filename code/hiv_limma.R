library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#Human HIV (3) dataset
#####

raw_counts_hiv <- "https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_ASVs_table.tsv"
conds_hiv <-"https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_metadata.tsv"

hiv <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
conditions_h <- read.table(file=conds_hiv, sep='\t', row.names = 1, header = T)

# filter to remove any ASV present in <10 % of samples,
hiv.filt <- hiv[apply(hiv, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

hiv.data.limma <- lim.fun(data = as.matrix(hiv.filt), conditions = conditions_h$comparison, nloop = 100, mean = 0, prop_null = 0.95)
save(hiv.data.limma, file="/Volumes/data2/andreea/ext_analysis/hiv.data.limma.Rda") 