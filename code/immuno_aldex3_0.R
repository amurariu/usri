# you will need to change the path to ALDEx3 repository
devtools::load_all('~/Documents/github/ALDEx3')
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/ald3.fun.R')

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
immuno.data_0.aldex3 <- ald3.fun(data=immuno.data, conds=immuno.conds$conditions_p, nloop=100, gamma=1e-3)
save(immuno.data_0.aldex3, file="/Volumes/data2/andreea/ext_analysis/immuno.data.aldex3_0.Rda")

# 
# # aldex3 test with selex
# url <- "https://raw.githubusercontent.com/ggloor/datasets/refs/heads/main/selex.txt"
# selex <- read.table(url, header=T, row.names=1)
# 
# condition <- c(rep(0, 7), rep(1, 7))
# X <- formula(~condition)
# data <- data.frame(condition=condition)
# 
# nsample <- 256 # use this nsamples
# foo0 <- aldex(selex, X, data=data, nsample=nsample, scale=clr.sm, gamma=1e-3)
# foo1 <- aldex(selex, X, data=data, nsample=nsample, scale=clr.sm, gamma=1)
# 
# # get the final stats
# # parameter - which pair is being testes
# # entity - the name of the gene or feature
# # estimate - LFC estimate
# # std.error - SE of the LFC (not the same as diff.win)
# # p.val.adj - adjusted p value
# sum.sel <- summary.aldex(foo0)
# sum.sel1 <- summary.aldex(foo1)
# # significant parts
# 
# tr0 <- which(sum.sel$p.val.adj < 0.05)
# tr1 <- which(sum.sel1$p.val.adj < 0.05)
# 
# plot(sum.sel$std.error*sqrt(14), sum.sel$estimate)
# points(sum.sel1$std.error*sqrt(14), sum.sel1$estimate, col='blue', pch=19, cex=0.5)
# 
# points(sum.sel$std.error[tr0]*sqrt(14), sum.sel$estimate[tr0], col='red')
# points(sum.sel$std.error[tr1]*sqrt(14), sum.sel$estimate[tr1], col='blue', pch=19, cex=0.5)
