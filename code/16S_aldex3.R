# processing 16S rRNA amplicon data and running differential abundance scale 
# benchmarking - ALDEx3

# Greg Gloor & Scott Dos Santos
# Last edited: 15th July 2026

# This code processes publicly available 16S rRNA amplicon sequencing data for
# thinning and randomisation. Data are taken from Bian et al. (2017). mSphere;
# 2 (5): e00327-17. This study looked at gut microbiome differences between a
# range of different age groups, from kindergarten-aged children to centenarians
# and compared OTUs present in the respective groups. This analysis will use a
# comparison of kindergarteners vs. middle-schoolers using ALDEx3.

#################################### setup ####################################

devtools::load_all("~/Documents/GitHub/ALDEx3/")
library(seqgendiff)

source("code/ald3.fun.R")
source("code/get_confusion.R")

# load input data frames of 16S rRNA amplicon data
rda <- paste0("https://github.com/amurariu/usri/raw/refs/heads/main/data/tianyi_",
              c("mid", "kin"), "_ft.Rda")

for(i in rda){
  tf <- tempfile(fileext = ".Rda")
  download.file(url = i, destfile = tf, mode = "wb")
  load(tf)
  unlink(tf)
  rm(tf)
}

# bind input data and make metadata objects
counts <- cbind(kin,mid)
ages <- data.frame(group = c(rep("K", 103), rep("M", 114)))

# filter OTUs with seqgendiff, choosing the 1000 most "median-expressed" genes
set.seed(123)
counts.sub <- select_counts(mat = as.matrix(counts), nsamp = ncol(counts), 
                            ngene = 400,  gselect = "max")

# check for and remove zero sum rows
any(rowSums(counts.sub) == 0) # returns FALSE
# counts.sub <- counts.sub[-which(rowSums(counts.sub) == 0),]
# any(rowSums(counts.sub) == 0) # returns FALSE

############################ thinning & con matrix ############################

# loop through gamma values from 0.001 to 0.5 as per other datasets and run the
# ALDEx3 function to generate 100 iterations of randomised & thinned data
gamma <- c(1e-3, seq(0.1, 0.5, 0.1))
gval <- seq(0, 5, 1)

for(g in 1:length(gamma)){
  
  # make object name based on gamma value
  name <- paste0("16S.data_", gval[g], ".aldex3")
 
  # run aldex3 function for randomising and thinning; std = 4 due to much
  # larger variance in the count data
  ald3.out <-  ald3.fun(data = counts.sub, conds = ages$group, nloop=100,
                        gamma = gamma[g], prop_null = 0.95, mean = 0,
                        std = 4, nsample = 256)
  
  # assign to new variable and save
  assign(x = name, value = ald3.out)
  save(list = name, file = paste0("~/Documents/GitHub/ext_analysis/16S.data.aldex3_", gval[g], ".Rda" ))
  
  # generate confusion matrix for the 100 iterations of thinned/randomised data
  # for the corresponding gamma value and save
  name.cm <- paste0("cm.16S.aldex3.", gval[g])
  assign(x = name.cm, value = get_confusion(ald3.out, prog = "ALDEx3", FDR = 0.05))
  save(list = name.cm, file = paste0("~/Documents/GitHub/usri/analysis/confusionMats/", name.cm, ".Rda" ))
  
}

################################# check plots #################################

# effect plot for non-randomised, non-thinned dataset (i.e. original data) at
# gamma = 0
plot(x = `16S.data_0.aldex3`$u.data$std.error * sqrt(ncol(counts.sub)),
     y = `16S.data_0.aldex3`$u.data$estimate, xlab = "Log dispersion within",
     ylab = "Log difference between")

sig.u <- which(`16S.data_0.aldex3`$u.data$p.val.adj <0.05)

points(x = `16S.data_0.aldex3`$u.data$std.error[sig.u] * sqrt(ncol(counts.sub)),
       y = `16S.data_0.aldex3`$u.data$estimate[sig.u], pch = 19, col = "red")

abline(a = 0, b = 1, lty = 2, col = "grey")
abline(a = 0, b = -1, lty = 2, col = "grey")
abline(h = 0, lty = 2, col = "blue")


# effect plot for randomised/thinned dataset at gamma = 0
plot(x = `16S.data_0.aldex3`$t.data[[1]]$std.error * sqrt(ncol(counts.sub)),
     y = `16S.data_0.aldex3`$t.data[[1]]$estimate, xlab = "Log dispersion within",
     ylab = "Log difference between")

sig.t <- which(`16S.data_0.aldex3`$t.data[[1]]$p.val.adj <0.05)

points(x = `16S.data_0.aldex3`$t.data[[1]]$std.error[sig.t] * sqrt(ncol(counts.sub)),
       y = `16S.data_0.aldex3`$t.data[[1]]$estimate[sig.t], pch = 19, col = "red")

abline(a = 0, b = 1, lty = 2, col = "grey")
abline(a = 0, b = -1, lty = 2, col = "grey")
abline(h = 0, lty = 2, col = "blue")


# effect plot for randomised/thinned dataset at gamma = 5
plot(x = `16S.data_5.aldex3`$t.data[[1]]$std.error * sqrt(ncol(counts.sub)),
     y = `16S.data_5.aldex3`$t.data[[1]]$estimate, xlab = "Log dispersion within",
     ylab = "Log difference between")

sig.t5 <- which(`16S.data_5.aldex3`$t.data[[1]]$p.val.adj <0.05)

points(x = `16S.data_5.aldex3`$t.data[[1]]$std.error[sig.t5] * sqrt(ncol(counts.sub)),
       y = `16S.data_5.aldex3`$t.data[[1]]$estimate[sig.t5], pch = 19, col = "red")

abline(a = 0, b = 1, lty = 2, col = "grey")
abline(a = 0, b = -1, lty = 2, col = "grey")
abline(h = 0, lty = 2, col = "blue")
