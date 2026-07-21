devtools::load_all('~/Documents/GitHub/ALDEx3')
library(seqgendiff, warn.conflicts=F)

source('code/ald3.fun.R')


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

# run aldex3 with gamma = 0.1, then save
pseudo.data_1.aldex3 <- ald3.fun(data = sc.pseudo, conds = conds.stim$stim, nloop = 100, gamma=0.1, prop_null = 0.95, mean = 0, std = 2)
save(pseudo.data_1.aldex3, file="~/Documents/GitHub/ext_analysis/pseudo_sc.data.aldex3_1.Rda")
