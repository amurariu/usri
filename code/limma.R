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

# immuno is the data table
# immuno.conds is the conditions for the unpermuted data
# N is the number of random instances
immuno.data.limma <- lim.fun(immuno.data, immuno.conds, 100)
save(immuno.data.limma, file="../ext_analysis/immuno.data.limma.Rda") 

brca.data.limma <- lim.fun(brca.data, brca.conds, 100)
save(brca.data.limma, file="../ext_analysis/brca.data.limma.Rda") 
