library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

# loads the limma.voom function to permute and add random true positives to the dataset
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

#####
#BRCA dataset
#####

#pull the data and filter
raw_counts_brca<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts_brca, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=conds_brca, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- as.vector(conditions_b) #changed from conditions to brca.conds for consistency with PD1 dataset

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
kirc.conds <- as.vector(unlist(read.table(file=conds_kirc, sep='\t'))) 

y_kirc <- DGEList(counts=kirc, group=factor(kirc.conds))
keep_kirc <- filterByExpr(y_kirc)
y_kirc <- y_kirc[keep_kirc,keep.lib.sizes=FALSE]
kirc.data <- y_kirc$counts #filtered base dataset

#####
#LIHC dataset
#####

raw_counts_lihc <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.normal-tumor.pair.rawCount.tsv"
conds_lihc <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LIHC.conditions.tsv"

lihc <- read.table(file=raw_counts_lihc, header=T, row.names=1, sep='\t')
lihc.conds <- as.vector(unlist(read.table(file=conds_lihc, sep='\t'))) 

y_lihc <- DGEList(counts=lihc, group=factor(lihc.conds))
keep_lihc <- filterByExpr(y_lihc)
y_lihc <- y_lihc[keep_lihc,keep.lib.sizes=FALSE]
lihc.data <- y_lihc$counts #filtered base dataset

#####
#LUAD dataset
#####

raw_counts_luad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.normal-tumor.pair.rawCount.tsv"
conds_luad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-LUAD.conditions.tsv"
luad <- read.table(file=raw_counts_luad, header=T, row.names=1, sep='\t')
luad.conds <- as.vector(unlist(read.table(file=conds_luad, sep='\t'))) 

y_luad <- DGEList(counts=luad, group=factor(luad.conds))
keep_luad <- filterByExpr(y_luad)
y_luad <- y_luad[keep_luad,keep.lib.sizes=FALSE]
luad.data <- y_luad$counts #filtered base dataset

#####
#PRAD dataset
#####

raw_counts_prad <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.normal-tumor.pair.rawCount.tsv"
conds_prad <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-PRAD.conditions.tsv"
prad <- read.table(file=raw_counts_prad, header=T, row.names=1, sep='\t')
prad.conds <- as.vector(unlist(read.table(file=conds_prad, sep='\t'))) 

y_prad <- DGEList(counts=prad, group=factor(prad.conds))
keep_prad <- filterByExpr(y_prad)
y_prad <- y_prad[keep_prad,keep.lib.sizes=FALSE]
prad.data <- y_prad$counts #filtered base dataset

#####
#THCA dataset
#####

raw_counts_thca <- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.normal-tumor.pair.rawCount.tsv"
conds_thca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-THCA.conditions.tsv"
thca <- read.table(file=raw_counts_thca, header=T, row.names=1, sep='\t')
thca.conds <- as.vector(unlist(read.table(file=conds_thca, sep='\t'))) 

y_thca <- DGEList(counts=thca, group=factor(thca.conds))
keep_thca <- filterByExpr(y_thca)
y_thca <- y_thca[keep_thca,keep.lib.sizes=FALSE]
thca.data <- y_thca$counts #filtered base dataset

# immuno is the data table
# immuno.conds is the conditions for the unpermuted data
# N is the number of random instances
immuno.data.limma <- lim.fun(immuno.data, immuno.conds, 100)
save(immuno.data.limma, file="../ext_analysis/immuno.data.limma.Rda") 

brca.data.limma <- lim.fun(brca.data, brca.conds, 100)
save(brca.data.limma, file="../ext_analysis/brca.data.limma.Rda") 

kirc.data.limma <- lim.fun(kirc.data, kirc.conds, 100)
save(kirc.data.limma, file="../ext_analysis/kirc.data.limma.Rda") 

lihc.data.limma <- lim.fun(lihc.data, lihc.conds, 100)
save(lihc.data.limma, file="../ext_analysis/lihc.data.limma.Rda") 

luad.data.limma <- lim.fun(luad.data, luad.conds, 100)
save(luad.data.limma, file="../ext_analysis/luad.data.limma.Rda") 

prad.data.limma <- lim.fun(prad.data, prad.conds, 100)
save(prad.data.limma, file="../ext_analysis/prad.data.limma.Rda") 

thca.data.limma <- lim.fun(thca.data, thca.conds, 100)
save(thca.data.limma, file="../ext_analysis/thca.data.limma.Rda") 

