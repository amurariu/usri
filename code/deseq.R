# code for generating DESeq2 data x100 (all datasets)

# Author: Andreea Murariu
# Last edited: 2025-06-16

library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

# loads the DESeq2 function to randomise groups and add random true positives to
# the count tables
source('code/des.fun.R')

#####
# PD1 immunotherapy dataset
####

# pull the raw count data from github
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'

immuno <- read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')

# make conditions vector for PD1 dataset using metadata info (Pre vs. On)
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p)

# edgeR conditions for initial filtering
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset

#####
#BRCA dataset
#####

#pull the data and filter as above
raw_counts<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
con <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=con, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- data.frame(conditions_b) #changed from conditions to brca.conds for consistency with PD1 dataset

y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset

# immuno is the data table
# immuno.conds is the vector of conditions for the non-randomised data
# nloop is the number of randomisation instances for randomised data
immuno.data.DESeq <- des.fun(data = immuno.data, nloop = 100, conditions = immuno.conds$conditions_p)

save(immuno.data.DESeq, file="../ext_analysis/immuno.data.deseq.Rda") 

#BRCA function
brca.data.DESeq <- des.fun(brca.data, brca.conds, 100)
save(brca.data.DESeq, file="../ext_analysis/brca.data.deseq.Rda")
