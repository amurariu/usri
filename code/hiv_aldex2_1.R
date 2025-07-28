library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

source('code/ald2.fun.R')

#####
#Human HIV (3) dataset
#####

raw_counts_hiv <- "https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_ASVs_table.tsv"
conds_hiv <-"https://raw.githubusercontent.com/amurariu/usri/main/data/hiv_noguerajulian_metadata.tsv"

hiv <- read.table(file=raw_counts_hiv, header=T, row.names=1, sep='\t')
conditions_h <- read.table(file=conds_hiv, sep='\t', row.names = 1, header = T)
hiv.conds <- data.frame(conditions_h) 

# filter to remove any ASV present in <10 % of samples,
hiv.filt <- hiv[apply(hiv, 1, function(x){length(which(x != 0))/length(x)} >0.1),]

hiv.data_1.aldex2 <- ald2.fun(data=as.matrix(hiv.filt), conditions=hiv.conds$comparison, nloop=100, gamma=0.1)
save(hiv.data_1.aldex2, file="/Volumes/data2/andreea/ext_analysis/hiv.data.aldex2_1.Rda")