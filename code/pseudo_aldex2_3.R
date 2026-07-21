library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

source('code/ald2.fun.R')


# load data from github
url.counts <- url("https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/pseudobulk_counts.txt")
sc.pseudo <- read.table(url.counts, header = T, sep = "\t", row.names = 1)

# remove zero sum rows
sc.pseudo <- sc.pseudo[-which(rowSums(sc.pseudo) == 0),]

# convert counts to matrix
sc.pseudo <- as.matrix(sc.pseudo)

# filter for top 8,500 median-expressed genes
set.seed(123)
sc.pseudo <- select_counts(mat = sc.pseudo, nsamp = ncol(sc.pseudo),
                           ngene = 8500, gselect = "max", filter_first = F)

# make conditions vector
conds.stim <- data.frame(stim = c(rep("CRTL", 103), rep("STIM", 103)))
rownames(conds.stim) <- colnames(sc.pseudo)

# run aldex2 with gamma = 0.3, then save
pseudo.data_3.aldex2 <- ald2.fun(data = sc.pseudo, conditions = conds.stim$stim, nloop=100, gamma=0.3)
save(pseudo.data_3.aldex2, file="~/Documents/GitHub/ext_analysis/pseudo_sc.data.aldex2_3.Rda")
