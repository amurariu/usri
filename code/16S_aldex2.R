# processing 16S rRNA amplicon data and running differential abundance scale 
# benchmarking - ALDEx2

# Greg Gloor & Scott Dos Santos
# Last edited: 6th July 2026

# This code processes publicly available 16S rRNA amplicon sequencing data for
# thinning and randomisation. Data are taken from Bian et al. (2017). mSphere;
# 2 (5): e00327-17. This study looked at gut microbiome differences between a
# range of different age groups, from kindergarten-aged children to centenarians
# and compared OTUs present in the respective groups. This analysis will use a
# comparison of kindergarteners vs. middle-schoolers.

#################################### setup ####################################

library(ALDEx2)
library(seqgendiff)

source("code/ald2.fun.R")
source("code/get_confusion.R")

# load and bind input data frames of 16S rRNA amplicon data
load("~/Documents/GitHub/tiyani-aldex3/data/kin.Rda") # private repo
load("~/Documents/GitHub/tiyani-aldex3/data/mid.Rda") # private repo

counts <- cbind(kin,mid)
ages <- data.frame(group = c(rep("K", 103), rep("M", 114)))

# filter OTUs with seqgendiff, choosing the 1000 most "median-expressed" genes
set.seed(123)
counts.sub <- select_counts(mat = as.matrix(counts), nsamp = ncol(counts), 
                            ngene = 1000,  gselect = "max")

# check for and remove zero sum rows
any(rowSums(counts.sub) == 0) # returns TRUE
counts.sub <- counts.sub[-which(rowSums(counts.sub) == 0),]
any(rowSums(counts.sub) == 0) # returns FALSE; removed 6 rows

############################ thinning & con matrix ############################

# loop through gamma values from 0.001 to 0.5 as per other datasets and run the
# ALDEx2 function to generate 100 iterations of randomised & thinned data
gamma <- c(1e-3, seq(0.1, 0.5, 0.1))
gval <- seq(0, 5, 1)

for(g in 1:length(gamma)){
  
  # make object name based on gamma value
  name <- paste0("16S.data_", gval[g], ".aldex2")
  
  # run aldex3 function for randomising and thinning
  ald2.out <-  ald2.fun(data = counts.sub, conditions = ages$group, nloop=100,
                        gamma = gamma[g], prop_null = 0.05, mean = 0,
                        std = 2, nsample = 128)
  
  # assign to new variable and save
  assign(x = name, value = ald2.out)
  save(list = name, file = paste0("~/Documents/GitHub/usri/data/16S.data.aldex2_", gval[g], ".Rda" ))
  
  # generate confusion matrix for the 100 iterations of thinned/randomised data
  # for the corresponding gamma value and save
  name.cm <- paste0("cm.16S.aldex2.", gval[g])
  assign(x = name.cm, value = get_confusion(ald2.out, prog = "ALDEx2", FDR = 0.05))
  save(list = name.cm, file = paste0("~/Documents/GitHub/usri/analysis/confusionMats/", name.cm, ".Rda" ))
  
}

################################# check plots #################################

# effect plot for non-randomised, non-thinned dataset (i.e. original data) at
# gamma = 0
plot(x = `16S.data_0.aldex2`$u.data$diff.win,
     y = `16S.data_0.aldex2`$u.data$diff.btw,
     xlab = "Log dispersion within", ylab = "Log difference between")

sig.u <- which(`16S.data_0.aldex2`$u.data$we.eBH <0.05)

points(x = `16S.data_0.aldex2`$u.data$diff.win[sig.u],
       y = `16S.data_0.aldex2`$u.data$diff.btw[sig.u], pch = 19, col = "red")

abline(a = 0, b = 1, lty = 2, col = "grey")
abline(a = 0, b = -1, lty = 2, col = "grey")
abline(h = 0, lty = 2, col = "blue")


# effect plot for randomised/thinned dataset at gamma = 0
plot(x = `16S.data_0.aldex2`$t.data[[1]]$diff.win,
     y = `16S.data_0.aldex2`$t.data[[1]]$diff.btw,
     xlab = "Log dispersion within", ylab = "Log difference between")

sig.t <- which(`16S.data_0.aldex2`$t.data[[1]]$we.eBH <0.05)

points(x = `16S.data_0.aldex2`$t.data[[1]]$diff.win[sig.t],
       y = `16S.data_0.aldex2`$t.data[[1]]$diff.btw[sig.t], pch = 19, col = "red")

abline(a = 0, b = 1, lty = 2, col = "grey")
abline(a = 0, b = -1, lty = 2, col = "grey")
abline(h = 0, lty = 2, col = "blue")


# effect plot for randomised/thinned dataset at gamma = 5
plot(x = `16S.data_5.aldex2`$t.data[[1]]$diff.win,
     y = `16S.data_5.aldex2`$t.data[[1]]$diff.btw,
     xlab = "Log dispersion within", ylab = "Log difference between")

sig.t5 <- which(`16S.data_5.aldex2`$t.data[[1]]$we.eBH <0.05)

points(x = `16S.data_5.aldex2`$t.data[[1]]$diff.win[sig.t5],
       y = `16S.data_5.aldex2`$t.data[[1]]$diff.btw[sig.t5], pch = 19, col = "red")

abline(a = 0, b = 1, lty = 2, col = "grey")
abline(a = 0, b = -1, lty = 2, col = "grey")
abline(h = 0, lty = 2, col = "blue")
