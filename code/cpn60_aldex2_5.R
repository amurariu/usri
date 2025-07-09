# running aldex2 (0.5 gamma) on cpn60 dataset (ST3, abx vs. no abx exposure)

# Scott Dos Santos
# 2025-06-23

library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

# load ALDEx2 function for generating randomised/thinned data
source('code/ald2.fun.R')

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

# convert count table to matrix
cpn60.counts <- as.matrix(cpn60.counts)

# generate results for thinned +/- randomised cpn60 data, 100 times, plus 1x
# non-randomised, non-thinned (i.e. original data)
cpn60.data.aldex2.0 <- ald2.fun(data = cpn60.counts, nloop = 100,
                                conditions = cpn60.meta$abx, gamma = 0.5)

save(cpn60.data.aldex2.0, file = "../ext_analysis/cpn60.data.aldex2_5.Rda")
