# running aldex2 (0.2 gamma) on vaginal metatranscriptome dataset (healthy vs. BV)

# Scott Dos Santos
# 2025-07-10

library(ALDEx2) #do not load aldex2 and aldex3 at the same time
library(seqgendiff)

# load ALDEx2 function for generating randomised/thinned data
source('code/ald2.fun.R')

# load metatranscriptome (mts) dataset and metadata
loc.counts <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/mts.data.filt.txt"
loc.meta <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/mts.metadata.txt"

mts.counts <- read.table(loc.counts, header = T, sep = "\t", quote = "", row.names = 1)
mts.meta <- read.table(loc.meta, header = T, sep = "\t", quote = "", row.names = 1)

# convert count table to matrix
mts.counts <- as.matrix(mts.counts)

# generate results for thinned +/- randomised mts data, 100 times, plus 1x
# non-randomised, non-thinned (i.e. original data)
mts.data.aldex2.2 <- ald2.fun(data = mts.counts, nloop = 100,
                              conditions = mts.meta$groups.2, gamma = 0.2)

save(mts.data.aldex2.2, file = "../ext_analysis/mts.data.aldex2_2.Rda")
