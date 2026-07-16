# processing single-cell data w/ pseudobulk apparoach for differential abundance
# scale benchmarking

# Greg Gloor & Scott Dos Santos
# Last edited: 30th June 2026

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
library(SeuratData)
library(Matrix)

# clone the seurat-data repo (https://github.com/satijalab/seurat-data) and then
# load via devtools
# devtools::load_all("~/Documents/GitHub/seurat-data/")

# check available seurat datasets
AvailableData()

# install the ifnb dataset; ONLY NEEDS TO BE RUN ONCE!
# InstallData('ifnb')

# load the ifnb dataset (ignore warnings)
ifnb.orig <- LoadData("ifnb")
ifnb <- ifnb.orig
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

##### added by GG July 16, 2026
# get metadata
meta <- ifnb@meta.data
# make a master vector of unique names STIM_DONOR_TYPE to use for aggregation
# substitute spaces with _ in seurat_annotation column
meta$seurat_annotations <- gsub(" ", "_", meta$seurat_annotations)
# remove SNG- from the Donor names
meta$donor_id <- gsub("SNG-", "", meta$donor_id)
# concatenate names
con_names <- paste(meta$stim, meta$donor_id, meta$seurat_annotations, sep="_")

counts <- as.data.frame(Matrix(ifnb$RNA$counts, sparse=F))

# finally aggregate 
ifnb.agg <- stats::aggregate(t(counts), by=list(con_names), FUN=sum)
rownames(ifnb.agg) <- ifnb.agg$Group.1
ifnb.agg$Group.1 <- NULL
####

# subset into  
# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, 
                                   group.by = c("stim", "donor_id", "seurat_annotations"))

# add cell type and stimulated condition to metadata
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
Idents(pseudo_ifnb) <- "celltype.stim"

# export normalised/aggregated count data to matrix
seu.agg.counts <- as.data.frame(Matrix(pseudo_ifnb$RNA$counts,sparse=F))

#### added by GG July 16, 2026
#  head(colSums(seu.agg.counts), n=10L)
#  head(colSums(ifnb.agg), n=10L)
# tldr: equivalent columns give exactly the same read count total and 
# the same count per gene
####

# move rownames to own column and move to front of df
counts$feature <- rownames(counts)
counts <- counts[,c(ncol(counts), 1:ncol(counts)-1)]

# write to txt
write.table(counts, file = paste0("~/Documents/GitHub/usri/data/pseudobulk_counts.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)

