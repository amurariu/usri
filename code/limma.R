library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

# loads the limma.voom function to permute and add random true positives to the dataset
source('code/lim.fun.R')

#####
# PD1 immunotherapy dataset
####

# pull the data and filter using edgeR
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p)

#edgeR conditions for initial filtering
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE] #had to take out y_pd1$counts
immuno.data <- y_pd1$counts

# immuno is the data table
# immuno.conds is the conditions for the unpermuted data
# N is the number of random instances
immuno.data.limma <- lim.fun(immuno.data, conditions_p, 10)
save(immuno.data.limma, file="../ext_analysis/immuno.data.limma.Rda") 
