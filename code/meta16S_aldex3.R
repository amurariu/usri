####
#
# Generates the aldex3 thinned data for the kin vs mid Tiyani dataset
# Greg Gloor
# 
# The trick is to use the built-in filtering from seqgendiff 
# via "select_counts", rather than the edgeR or codaSeq filtering
#
###

load('~/Documents/0_git/projects/tiyani-aldex3/data/kin.Rda')
load('~/Documents/0_git/projects/tiyani-aldex3/data/mid.Rda')

km <- cbind(kin, mid)

library(ALDEx3)
library(seqgendiff, warn.conflicts=F)

source('~/Documents/0_git/projects/andreea_usri/code/ald3.fun.R')
source('~/Documents/0_git/projects/andreea_usri/code/get_confusion.R')

conds.km <- c(rep("K", 103), rep("M", 114))

# seqgendiff code to subset matrix for thinning
# keeping 1000 of 1117 OTUs
sub.med <- select_counts(as.matrix(km), nsamp=217, 1000, gselect='max')

# keeps all 48 features that have > 0
# none of the remainder are only 0
# which(apply(jnk.med, 1, min) > 0)
# length(which(apply(jnk.med, 1, sum) < 100))
# returns 0

# thinning works on this subset!!

gamma = c(1e-3, 0.1,0.2,0.3,0.4,0.5)
gname = c(0,1,2,3,4,5)
location <- "~/Documents/0_git/projects/andreea_usri/"
for(g in 1:length(gamma)){

nm.ald3 <- paste("meta16S_gamma",gname[g],"aldex3", sep="")

# using smaller prop null because microbiome data is so variable
ald3.out <-  ald3.fun(data=sub.med, conds=conds.km, nloop=100, gamma=gamma[g], prop_null=0.05, mean=0, std=2, nsample=256)
assign(nm.ald3, ald3.out)

save(list=nm.ald3, file=paste(location, "data/", nm.ald3, ".Rda", sep="") )

# this makes a character
nm.cm <- paste0("cm.meta16S.aldex3.",gname[g], sep="")

#this makes a variable with the name in nm.cm
assign(nm.cm, get_confusion(ald3.out, prog="ALDEx3"))
save(list=nm.cm, file=paste(location, "analysis/confusionMats/", nm.cm, ".Rda", sep="") )
}

load("data/meta16S_gamma0aldex3.Rda")
load("analysis/confusionMats//cm.meta16S.aldex3.0.Rda")

cm.meta16S.aldex3.0$TPFPR

# original dataset
plot(meta16S_gamma0aldex3$u.data$std.error*sqrt(218), meta16S_gamma0aldex3$u.data$estimate)

# randomized/thinned dataset
plot(meta16S_gamma0aldex3$t.data[[1]]$std.error*sqrt(218), meta16S_gamma0aldex3$t.data[[1]]$estimate)
  