###
# June 25, 2026 by gg
#
# This code generates a pseudo-bulk RNA-seq dataset from single cell dataset
# found as an exemplar in the Seurat R package for ss RNA-seq data
# not surprisingly, it turns out the a lot of ssRNA-seq stuff is built off of
# repackaged bulk RNA seq and then hidden away in an aggregator package
# This makes it easy for people to use the tools but makes it hard
# to know exactly what is going on unless you dig a bit
# 
#####
# NOTE: the stuff you need is all at the top, how I got there is commented out
# and below the #### getting data #### line
### 
## https://satijalab.org/seurat/articles/de_vignette
###
## how to get Seurat data
## install Seurat from CRAN
## install.packages("Seurat")
## install seurat_data
##   maddeningly, this is not in CRAN or BioConductor but is available 
##   at https://github.com/satijalab/seurat-data
##   equally maddeningly, the install instructions did not work for me
##   so I cloned the repository to my local
###

###
# data are from Multiplexed droplet single-cell RNA-sequencing using natural genetic variation
# Nature Biotech 36:89,  2018
# referred to as the ifnb dataset in Seurat DE vignette
###

library(ggplot2)
library(Seurat)
# devtools::load_all('~/path/to/seurat-data')
 devtools::load_all('~/Documents/0_git/projects/seurat-data')
library(Matrix)
devtools::load_all('~/Documents/0_git/projects/ALDEx3')
library(ALDEx2)

# loads pre-computed ALDEx2 output
load(url('https://github.com/amurariu/usri/blob/main/data/ss.bulk_x.all.Rda?raw=TRUE') )
# ss.bulk.df
load(url('https://github.com/amurariu/usri/blob/main/data/ss.bulk.Rda?raw=TRUE'))
####
# ALDEx2
####
# extract out the pseudo-bulk data
# x <- aldex.clr(ss.bulk.df, cond=cond.stim)
# x.e <- aldex.effect(x)
# x.t <- aldex.ttest(x)
# ss.bulk_x.all <- cbind(x.e,x.t)
# 
# save(ss.bulk_x.all, file='~/Documents/0_git/projects/andreea_usri/data/ss.bulk_x.all.Rda') 

######
# plot(x.all) shows the stimulated cells have many expressed genes not in
# non-stimulated. Essentially nothing is down-regulated
# great example of asymmetric expression
######

####
# ALDEx3 total pseudobulk
####
conds.stim <- c(rep("CRTL", 103), rep("STIM", 103))
Y <- as.matrix(ss.bulk.df)
# Metadata Samples x Metadata Values
X <- data.frame(conds.stim)
X$conds.stim <- factor(X$conds.stim)

# Relative Abundances CRTL vs. STIM
# need small sample numbers for plotting effect sizes
#aldex.tss <- ALDEx3::aldex(Y, ~conds.stim, X, nsample=16, scale=tss.sm, gamma=1e-3) 
aldex.ss <- ALDEx3::aldex(Y, ~conds.stim, X, nsample=16, scale=clr.sm, gamma=0.3) 

ALDEx3::aldex.plot(aldex.ss, contrast='conds.stim', plot='effect', cex=0.5)
abline(h=0, lty=2, col='cyan', lwd=1.5)

###########
## pseudobulk Naive CD4
###########
naive.conds <- c(rep("CRTL", 8), rep("STIM", 8))
naive4 <- Y[,grep("Naive", colnames(Y))]
X4 <- data.frame(naive.conds)
X4$naive.conds <- factor(X4$naive.conds)
aldex.sb <- ALDEx3::aldex(naive4, ~naive.conds, X4, nsample=128, scale=clr.sm, gamma=1e-3)
ALDEx3::aldex.plot(aldex.sb, contrast='naive.conds', plot='effect')
abline(h=0, lty=2, col='cyan', lwd=1.5)



########### getting data ##########
# # commented out, but uncomment first character to get runnable code
# ## you might think that the vignette on DA would work
# ## wrong!!!, this is massaged and cleaned up
# 
# AvailableData() # shows available Seurat datasets
# # InstallData('ifnb') # installs the ifnb dataset-go for a walk
# 
# ## so now the vignette will work?
# ifnb <- LoadData("ifnb")
# 
# # Normalize the data
# # both work 
# # ifnb <- NormalizeData(ifnb) # this is proportion scaled by 10000
# ifnb <- NormalizeData(ifnb, normalization.method='CLR')
# 
# # Find DE features between CD16 Mono and CD14 Mono
# Idents(ifnb) <- "seurat_annotations"
# monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono")
# # view results
# head(monocyte.de.markers)
# 
# ### pseudobulk
# # load the inferred sample IDs of each cell
# ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
# stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
# info <- rbind(ctrl, stim)
# 
# # rename the cell IDs by substituting the '-' into '.'
# info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)
# 
# # only keep the cells with high-confidence sample ID
# info <- info[grep(pattern = "SNG", x = info$BEST), ]
# 
# # remove cells with duplicated IDs in both ctrl and stim groups
# info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]
# 
# # now add the sample IDs to ifnb 
# rownames(info) <- info$BARCODE
# info <- info[, c("BEST"), drop = F]
# names(info) <- c("donor_id")
# ifnb <- AddMetaData(ifnb, metadata = info)
# 
# # remove cells without donor IDs
# ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
# ifnb <- subset(ifnb, subset = donor_id != "unknown")
# 
# # pseudobulk the counts based on donor-condition-celltype
# pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
# 
# pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
# 
# Idents(pseudo_ifnb) <- "celltype.stim"
# 
# bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
#                          ident.1 = "CD14 Mono_STIM", 
#                          ident.2 = "CD14 Mono_CTRL",
#                          test.use = "DESeq2")
# head(bulk.mono.de, n = 15)
# 
# ss.bulk <- Matrix(pseudo_ifnb$RNA$counts,sparse=F)
# ss.bulk.df <- as.data.frame(ss.bulk)
# conds.stim <- c(rep("CRTL", 103), rep("STIM", 103))
# save(ss.bulk.df, file='~/Documents/0_git/projects/andreea_usri/data/ss.bulk.Rda')
# 
