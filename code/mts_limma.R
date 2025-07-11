# running Limma on vaginal metatranscriptome dataset (healthy vs. BV)

# Scott Dos Santos
# 2025-07-11

library(limma)
library(seqgendiff)

# load Limma function for generating randomised/thinned data
source('code/lim.fun.R')

# load metatranscriptome (mts) dataset and metadata
loc.counts <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/mts.data.filt.txt"
loc.meta <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/mts.metadata.txt"

mts.counts <- read.table(loc.counts, header = T, sep = "\t", quote = "", row.names = 1)
mts.meta <- read.table(loc.meta, header = T, sep = "\t", quote = "", row.names = 1)

# convert count table to matrix
mts.counts <- as.matrix(mts.counts)

# generate results for thinned +/- randomised mts data, 100 times, plus 1x
# non-randomised, non-thinned (i.e. original data)
mts.data.limma <- lim.fun(data = mts.counts, conditions = mts.meta$groups.2,
                          nloop = 2, mean = 0, prop_null = 0.95)

save(mts.data.limma, file = "../ext_analysis/mts.data.limma.Rda")
