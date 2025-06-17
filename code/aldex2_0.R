library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald2.fun.R')

#immuno/PD1 dataset loading
raw_counts_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts_immuno, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta_immuno, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p) #changed conditions to conditions_p to be consistent across datasets

#edgeR conditions for initial filtering
#PD1
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset

#save file - gamma=1e-3
#immuno.data_0.aldex2 <- ald2.fun(data=immuno.data, conditions=immuno.conds, nloop=100, gamma=1e-3)
#save(immuno.data_0.aldex2, file="../ext_analysis/immuno.data.aldex2_0.out.Rda")

#save file - gamma=0.2
#immuno.data_2.aldex2 <- ald2.fun(data=immuno.data, conditions=immuno.conds, nloop=100, gamma = 0.2)
#save(immuno.data_2.aldex2, file="../ext_analysis/immuno.data.aldex2_2.out.Rda")

#save file - gamma=0.5
#immuno.data_5.aldex2 <- ald2.fun(data=immuno.data, conditions=immuno.conds, nloop=100, gamma = 0.5)
#save(immuno.data_5.aldex2, file="../ext_analysis/immuno.data.aldex2_5.out.Rda")

#####
#BRCA dataset
#####

#pull the data and filter
raw_counts_brca<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
conds_brca <-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts_brca, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=conds_brca, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- data.frame(conditions_b) #changed from conditions to brca.conds for consistency with PD1 dataset

#edgeR
y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset

#save file - gamma=1e-3
brca.data_0.aldex2 <- ald2.fun(data=brca.data, conditions=brca.conds, nloop=100, gamma=1e-3)
save(brca.data_0.aldex2, file="../ext_analysis/brca.data.aldex2_0.out.Rda")

#save file - gamma=0.2
brca.data_2.aldex2 <- ald2.fun(data=brca.data, conditions=brca.conds, nloop=100, gamma = 0.2)
save(brca.data_2.aldex2, file="../ext_analysis/brca.data.aldex2_2.out.Rda")

#save file - gamma=0.5
brca.data_5.aldex2 <- ald2.fun(data=brca.data, conditions=brca.conds, nloop=100, gamma = 0.5)
save(brca.data_5.aldex2, file="../ext_analysis/brca.data.aldex2_5.out.Rda")

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

#save file - gamma=1e-3
kirc.data_0.aldex2 <- ald2.fun(data=kirc.data, conditions=kirc.conds, nloop=100, gamma=1e-3)
save(kirc.data_0.aldex2, file="../ext_analysis/kirc.data.aldex2_0.out.Rda")

#save file - gamma=0.2
kirc.data_2.aldex2 <- ald2.fun(data=kirc.data, conditions=kirc.conds, nloop=100, gamma = 0.2)
save(kirc.data_2.aldex2, file="../ext_analysis/kirc.data.aldex2_2.out.Rda")

#save file - gamma=0.5
kirc.data_5.aldex2 <- ald2.fun(data=kirc.data, conditions=kirc.conds, nloop=100, gamma = 0.5)
save(kirc.data_5.aldex2, file="../ext_analysis/kirc.data.aldex2_5.out.Rda")


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

#save file - gamma=1e-3
lihc.data_0.aldex2 <- ald2.fun(data=lihc.data, conditions=lihc.conds, nloop=100, gamma=1e-3)
save(lihc.data_0.aldex2, file="../ext_analysis/lihc.data.aldex2_0.out.Rda")

#save file - gamma=0.2
lihc.data_2.aldex2 <- ald2.fun(data=lihc.data, conditions=lihc.conds, nloop=100, gamma = 0.2)
save(lihc.data_2.aldex2, file="../ext_analysis/lihc.data.aldex2_2.out.Rda")

#save file - gamma=0.5
lihc.data_5.aldex2 <- ald2.fun(data=lihc.data, conditions=lihc.conds, nloop=100, gamma = 0.5)
save(lihc.data_5.aldex2, file="../ext_analysis/lihc.data.aldex2_5.out.Rda")
