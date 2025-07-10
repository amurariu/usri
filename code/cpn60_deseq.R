# running DESeq2 on cpn60 dataset (ST3, abx vs. no abx exposure)

# Scott Dos Santos
# 2025-06-25

library(seqgendiff, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

# load DESeq2 function for generating randomised/thinned data
source('code/des.fun.R')

# load cpn60 dataset and metadata
loc.counts <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/cpn60.data.txt"
loc.meta <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/cpn60.metadata.txt"

cpn60.counts <- read.table(loc.counts, header = T, sep = "\t", quote = "", row.names = 1)
cpn60.meta <- read.table(loc.meta, header = T, sep = "\t", quote = "", row.names = 1)

# substitute periods for dashes in column names of table
colnames(cpn60.counts) <- gsub("\\.", "-", colnames(cpn60.counts))

# filter for only ST3 samples (10 day-old infant stool)
cpn60.meta <- cpn60.meta[which(cpn60.meta$age == "10 days"),]
cpn60.counts <- cpn60.counts[, rownames(cpn60.meta)]

# check for and remove ASVs with zero sums
any(rowSums(cpn60.counts) == 0) # returns TRUE
cpn60.counts <- cpn60.counts[-which(rowSums(cpn60.counts)==0),]

# convert count table to matrix and add 1 to all counts (DESeq throws an error
# if every feature has a count of zero in one or more samples). This not at all
# statistically correct, but is nevertheless common practice (for example, see:
# https://www.biostars.org/p/440379/#9559966 
# https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564
# https://www.reddit.com/r/bioinformatics/comments/vxxgm2/having_an_error_for_days/
# https://github.com/cafferychen777/ggpicrust2/issues/130 )
cpn60.counts <- as.matrix(cpn60.counts)
cpn60.counts <- cpn60.counts +1

# generate results for thinned +/- randomised cpn60 data, 100 times, plus 1x
# non-randomised, non-thinned (i.e. original data)
cpn60.data.DESeq <- des.fun(data = cpn60.counts, nloop = 100, conditions = cpn60.meta$abx)
save(cpn60.data.DESeq, file = "../ext_analysis/cpn60.data.deseq.Rda")
