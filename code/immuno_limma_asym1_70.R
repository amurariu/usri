library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
# PD1 immunotherapy dataset
####

# pull the data and filter using edgeR
raw_counts_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts_immuno, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta_immuno, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109) #use str to figure out what conds p is
conditions_p[grep("_On",m)] <- "On"
immuno.conds<- as.vector(conditions_p)

#edgeR conditions for initial filtering
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE] #had to take out y_pd1$counts
immuno.data <- y_pd1$counts

immuno.data.limma.asym1.70 <- lim.fun(immuno.data, immuno.conds, 100, mean = 1, prop_null = 0.70)
save(immuno.data.limma.asym1.70, file="/Volumes/data2/andreea/ext_analysis/immuno.data.limma.asym1.prop70.Rda") 
