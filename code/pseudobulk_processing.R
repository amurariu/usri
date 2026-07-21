# processing single-cell data w/ pseudobulk apparoach for differential abundance
# scale benchmarking

# Greg Gloor & Scott Dos Santos
# Last edited: 21st July 2026

# This code generates a pseudo-bulk RNA-seq dataset from single cell dataset
# found as an exemplar in the Seurat R package for scRNA-seq data. This dataset
# is taken from Kang et al. (2018). Nature Biotech. 36; 89-94. Data represent a
# single-cell RNA-seq experiment where peripheral blood mononuclear cells were
# stimulated or not with interferon beta and sequenced. The end result is 103
# pseudobulked samples in each group.

# NOTE: installing the 'seurat-data' R package was difficult; the package did
#       not seem to be available on bioconductor or on CRAN, so had to clone it
#       from github (https://github.com/satijalab/seurat-data) as the repo's
#       instructions for installing it via devtools did not work

#################################### setup ####################################

library(Seurat)
library(Matrix)
library(SeuratData) # if loading from cloned github repo, ignore & use devtools

# clone the seurat-data repo (https://github.com/satijalab/seurat-data) and then
# load via devtools
# devtools::load_all("~/Documents/GitHub/seurat-data/")

# check available seurat datasets and install if missing, otherwise, load the
# ifnb dataset (ignore the warnings!)
if(AvailableData()["ifnb.SeuratData","Installed"] == FALSE){
  AvailableData()
  
  InstallData('ifnb')
} else{
  
  ifnb.orig <- LoadData("ifnb")
  ifnb <- ifnb.orig
  
}

# load the inferred sample IDs of each cell from the Kang dataset github repo
# (origin: https://github.com/yelabucsf/demuxlet_paper_code/tree/master/fig3 but
# these files have been saved in their original form in the current study repo)
url.id.ctrl <- url("https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/pseudobulk_ye1.ctrl.8.10.sm.best")
url.id.stim <- url("https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/pseudobulk_ye2.stim.8.10.sm.best")

ctrl <- read.table(url.id.ctrl, head = T, stringsAsFactors = F)
stim <- read.table(url.id.stim, head = T, stringsAsFactors = F)

# bind dfs for control and stimulated cells
info <- rbind(ctrl, stim)

# rename the cell IDs by substituting the '-' into '.' so they match the ifnb
# object
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)

# only keep the cells with high-confidence sample ID
info <- info[grep(pattern = "SNG", x = info$BEST), ]

# remove cells with duplicated IDs in both ctrl and stim groups
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]

############################ normalise & aggregate ############################

# normalise data using a centred log-ratio
ifnb <- NormalizeData(ifnb, normalization.method='CLR')

# add the unique cell IDs to the ifnb seurat object
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
ifnb <- AddMetaData(ifnb, metadata = info)

# remove cells without donor IDs
ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
ifnb <- subset(ifnb, subset = donor_id != "unknown")

# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, 
                                   group.by = c("stim", "donor_id", "seurat_annotations"))

# add cell type and stimulated condition to metadata
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
Idents(pseudo_ifnb) <- "celltype.stim"

# export normalised/aggregated count data to matrix
seu.agg.counts <- as.data.frame(Matrix(pseudo_ifnb$RNA$counts,sparse=F))

# move rownames to own column and move to front of df
seu.agg.counts$feature <- rownames(seu.agg.counts)
seu.agg.counts <- seu.agg.counts[,c(ncol(seu.agg.counts), 1:ncol(seu.agg.counts)-1)]

# replace spaces with underscores in column names
colnames(seu.agg.counts) <- gsub(" ", "_", colnames(seu.agg.counts))

# write to txt
write.table(counts, file = paste0("~/Documents/GitHub/usri/data/pseudobulk_counts.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)

############################## check aggregation ##############################

# want to check that seurat's aggregation function is *ACTUALLY* returning raw
# counts instead of some sort of normalised counts. Manually aggregate counts

# re-process raw data before CLR normalisation
ifnb.manual <- ifnb.orig
ifnb.manual <- AddMetaData(ifnb.manual, metadata = info)
ifnb.manual$donor_id[is.na(ifnb.manual$donor_id)] <- "unknown"
ifnb.manual <- subset(ifnb.manual, subset = donor_id != "unknown")

# get metadata from seurat object
meta <- ifnb.manual@meta.data

# make a master vector of unique names following pattern of 'GROUP_DONOR_TYPE'
# to use for aggregation: start by substituting spaces with underscores in the
# 'seurat_annotation' column
meta$seurat_annotations <- gsub(" ", "_", meta$seurat_annotations)

# add back in the dash betwen 'SNG' and donor id, then concatenate sample names
meta$donor_id <- gsub("SNG_", "SNG-", meta$donor_id)
con_names <- paste(meta$stim, meta$donor_id, meta$seurat_annotations, sep="_")

# export as matrix and aggregate raw counts manually (takes a little while)
counts <- as.data.frame(Matrix(ifnb$RNA$counts, sparse=F))
ifnb.manual.agg <- stats::aggregate(t(counts), by=list(con_names), FUN=sum)

# move aggregated sample names from their own column to row names
rownames(ifnb.manual.agg) <- ifnb.manual.agg$Group.1
ifnb.manual.agg$Group.1 <- NULL

# transpose to get samples as columns and re-order based on seurat-aggregated
# data frame column names (excluding 'feature column')
ifnb.manual.agg <- as.data.frame(t(ifnb.manual.agg))
ifnb.manual.agg <- ifnb.manual.agg[,colnames(seu.agg.counts[2:207])]

# are the sums of each sample identical?
all(colSums(ifnb.manual.agg)==colSums(seu.agg.counts[,2:ncol(seu.agg.counts)])) # returns TRUE

# are the data frames identical?
all(ifnb.manual.agg == seu.agg.counts[,2:ncol(seu.agg.counts)]) # returns TRUE

# samples give exactly the same read counts per sample AND per gene
