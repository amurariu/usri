library(CoDaSeq)

# load 16S data from HMP oral microbiomes (keratinised gingiva vs. oral plaque;
# included with CoDaSeq package)
oral.data <- ak_op
oral.taxa <- hmpgenera

# add "otu." to rownames of both dfs
rownames(oral.data) <- paste0("otu.", rownames(oral.data))
rownames(oral.taxa) <- paste0("otu.", rownames(oral.taxa))

# filter for OTU present in >10 % samples
oral.filt <- oral.data[apply(oral.data, MARGIN = 1, function(x){(length(which(x != 0))/length(x))>0.1}),]

# filter OTU tax table for only those present in feature table
oral.taxa <- oral.taxa[rownames(oral.filt),]

# make df of conditions and total/filtered read counts and their differences
oral.meta <- data.frame(group = gsub("_.*", "", colnames(oral.data)),
                        reads.all = colSums(oral.data),
                        reads.filt = colSums(oral.filt))
oral.meta$reads.diff <- oral.meta$reads.all - oral.meta$reads.filt

oral.meta$group <- gsub("ak", "Keratinised gingiva", oral.meta$group)
oral.meta$group <- gsub("op", "Oral plaque", oral.meta$group)

# write dfs to files
write.table(oral.filt, file = "~/Documents/PostDoc_Western/Data/diffAbund/datasets/16S/oral_data.txt",
            sep = "\t", quote = F, row.names = T, col.names = T)
write.table(oral.taxa, file = "~/Documents/PostDoc_Western/Data/diffAbund/datasets/16S/oral_taxa.txt",
            sep = "\t", quote = F, row.names = T, col.names = T)
write.table(oral.meta, file = "~/Documents/PostDoc_Western/Data/diffAbund/datasets/16S/oral_meta.txt",
            sep = "\t", quote = F, row.names = T, col.names = T)