library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

source('code/edg.fun.R')

##PD1 immuno dataset

#loading and filtering of data
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p) 

#edgeR conditions for initial filtering
#PD1
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset


#####
#BRCA dataset
#####

#pull the data and filter
raw_counts<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
con <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=con, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- data.frame(conditions_b) #changed from conditions to brca.conds for consistency with PD1 dataset

#edgeR
y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset


#saving file
immuno.data.edgeR <-edg.fun(immuno.data, conditions_p,100)
save(immuno.data.edgeR, file="../ext_analysis/immuno.data.edger.Rda")

brca.data.edgeR <- edg.fun(brca.data, brca.conds, 100)
save(brca.data.edgeR, file="../ext_analysis/brca.data.edger.Rda")
