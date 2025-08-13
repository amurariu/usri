library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)

source('code/lim.fun.R')

#####
#yeast transcriptome
#####

##setting up data
url <- "https://raw.githubusercontent.com/amurariu/usri/main/data/transcriptome.tsv"
yeast <- read.table(url, header=T, row.names=1)

# remove the one gene with 0 reads
yeast <- yeast[rownames(yst) != "YOR072W-B",]

# Gierlinski:2015aa
yeast[,c('SNF2.6', 'SNF2.13','SNF2.25','SNF2.35')] <- NULL
yeast[,c('WT.21','WT.22','WT.25','WT.28','WT.34','WT.36')] <- NULL

##setting up conditions vector
conditions_y <- c(rep('S', 44), rep('W', 42))
yeast.conds <- data.frame(conditions_y)

# edgeR conditions for initial filtering
y_yeast <- DGEList(counts=yeast, group=factor(conditions_y))
keep_yeast <- filterByExpr(y_yeast)
y_yeast <- y_yeast[keep_yeast,keep.lib.sizes=FALSE]
yeast.data <- y_yeast$counts #filtered base dataset

yeast.data.limma <- lim.fun(data = yeast.data, conditions = yeast.conds$conditions_y, nloop = 100, mean = 0, prop_null = 0.95)
save(yeast.data.limma, file="/Volumes/data2/andreea/ext_analysis/yeast.data.limma.Rda") 
