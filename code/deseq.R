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
raw_counts_brca<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=con, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- data.frame(conditions_b) 

y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset

#####
#KIRC
#####

raw_counts_kirc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.normal-tumor.pair.rawCount.tsv"
conds_kirc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-KIRC.conditions.tsv"

kirc <- read.table(file=raw_counts_kirc, header=T, row.names=1, sep='\t')
conditions_k <- as.vector(unlist(read.table(file=conds_kirc, sep='\t'))) 
kirc.conds <- data.frame(conditions_k)

y_kirc <- DGEList(counts=kirc, group=factor(kirc.conds))
keep_kirc <- filterByExpr(y_kirc)
y_kirc <- y_kirc[keep_kirc,keep.lib.sizes=FALSE]
kirc.data <- y_kirc$counts

#####
#LIHC dataset
#####

raw_counts_lihc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.normal-tumor.pair.rawCount.tsv"
conds_lihc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.conditions.tsv"
lihc <- read.table(file=raw_counts_lihc, header=T, row.names=1, sep='\t')
conditions_li <- as.vector(unlist(read.table(file=conds_lihc, sep='\t'))) 
lihc.conds <- data.frame(conditions_li) 

y_lihc <- DGEList(counts=lihc, group=factor(lihc.conds))
keep_lihc <- filterByExpr(y_lihc)
y_lihc <- y_lihc[keep_lihc,keep.lib.sizes=FALSE]
lihc.data <- y_lihc$counts 

#####
#LUAD dataset
#####

raw_counts_luad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.normal-tumor.pair.rawCount.tsv"
conds_luad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.conditions.tsv"
luad <- read.table(file=raw_counts_luad, header=T, row.names=1, sep='\t')
conditions_lu <- as.vector(unlist(read.table(file=conds_luad, sep='\t'))) 
luad.conds <- data.frame(conditions_lu) 

y_luad <- DGEList(counts=luad, group=factor(luad.conds))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts 

#####
#PRAD dataset
#####

raw_counts_prad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.normal-tumor.pair.rawCount.tsv"
conds_prad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.conditions.tsv"
prad <- read.table(file=raw_counts_prad, header=T, row.names=1, sep='\t')
conditions_pr <- as.vector(unlist(read.table(file=conds_prad, sep='\t')))
prad.conds <- data.frame(conditions_pr) 

y_prad <- DGEList(counts=prad, group=factor(prad.conds))
keep_prad <- filterByExpr(y_prad)
y_prad <- y_prad[keep_prad,keep.lib.sizes=FALSE]
prad.data <- y_prad$counts 

#####
#THCA dataset
#####

raw_counts_thca <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.normal-tumor.pair.rawCount.tsv"
conds_thca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.conditions.tsv"
thca <- read.table(file=raw_counts_thca, header=T, row.names=1, sep='\t')
conditions_t <- as.vector(unlist(read.table(file=conds_thca, sep='\t')))
thca.conds <- data.frame(conditions_t) 

y_thca <- DGEList(counts=thca, group=factor(thca.conds))
keep_thca <- filterByExpr(y_thca)
y_thca <- y_thca[keep_thca,keep.lib.sizes=FALSE]
thca.data <- y_thca$counts 

# immuno is the data table
# immuno.conds is the vector of conditions for the non-randomised data
# nloop is the number of randomisation instances for randomised data

#Immuno function
immuno.data.DESeq <- des.fun(data = immuno.data, nloop = 100, conditions = immuno.conds$conditions_p)
save(immuno.data.DESeq, file="../ext_analysis/immuno.data.deseq.Rda") 

#BRCA function
brca.data.DESeq <- des.fun(data = brca.data, nloop = 100, conditions = brca.conds$conditions_b) 
save(brca.data.DESeq, file="../ext_analysis/brca.data.deseq.Rda")

#KIRC function
kirc.data.DESeq <- des.fun(data = kirc.data, nloop = 100, conditions = kirc.conds$conditions_k)
save(kirc.data.DESeq, file="../ext_analysis/kirc.data.deseq.Rda") 

#LIHC function
lihc.data.DESeq <- des.fun(data = lihc.data, nloop = 100, conditions = lihc.conds$conditions_li)
save(lihc.data.DESeq, file="../ext_analysis/lihc.data.deseq.Rda")

#LUAD function
luad.data.DESeq <- des.fun(data = luad.data, nloop = 100, conditions = luad.conds$conditions_lu) 
save(luad.data.DESeq, file="../ext_analysis/luad.data.deseq.Rda")

#PRAD function
prad.data.DESeq <- des.fun(data = prad.data, nloop = 100, conditions = prad.conds$conditions_pr)
save(prad.data.DESeq, file="../ext_analysis/prad.data.deseq.Rda") 

#THCA function
thca.data.DESeq <- des.fun(data = thca.data, nloop = 100, conditions = thca.conds$conditions_t)
save(thca.data.DESeq, file="../ext_analysis/thca.data.deseq.Rda")