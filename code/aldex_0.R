library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald.fun.R')

#immuno/PD1 dataset loading
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')
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

#save file
immuno.data.aldex2_0.out <- ald.fun(immuno.data, conditions_p, 2, gamma=1e-3)
save(immuno.data_0.aldex2, file="../ext_analysis/immuno.data.aldex2_0.out.Rda")

#save file
immuno.data.aldex2_0.out <- ald.fun(immuno.data, conditions_p, 2, gamma=0.2)
save(immuno.data_2.aldex2, file="../ext_analysis/immuno.data.aldex2_2.out.Rda")

#save file
immuno.data.aldex2_0.out <- ald.fun(immuno.data, conditions_p, 2, gamma=0.5)
save(immuno.data_5.aldex2, file="../ext_analysis/immuno.data.aldex2_5.out.Rda")

