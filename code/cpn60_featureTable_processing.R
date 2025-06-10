# feature table processing: cpn60 amplicon - stool samples from LEGACY

# Scott Dos Santos
# 2025-06-09

# takes the raw feature table and relevant lookup tables for a cpn60 amplicon
# sequencing dataset for 1030 stool samples and filters out low-abundance and
# prevalence ASVs (see: https://doi.org/10.3389/fcimb.2023.1144254 for paper
# and https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP368054 for SRA access) 

#################################### setup ####################################

library(dplyr)
library(CoDaSeq)
library(stringr)

# links for pulling files from github
loc.ft <- "https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/cpn60_featureTable.txt"
loc.tax <- "https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/cpn60_taxonomy.txt"

# load in feature table of read counts and corresponding taxonomy table
# NOTE: taxonomy table has been edited to remove ASVs which: 1) are eukaryotic;
#       2) have a confidence score of <0.8 at the phylum rank (non-target)
cpn.ft <- read.table(loc.ft, header = T, sep = "\t", quote = "", row.names = 1)
cpn.tax <- read.table(loc.tax, header = T, sep = "\t", quote = "", row.names = 1)

################################# filtering FT #################################

# remove the eukaryotic and non-target ASVs from feature table
cpn.ft.filt <- cpn.ft[rownames(cpn.tax),]

# the original study highlighted a single ASV likely representing a contaminant,
# Pseudomonas tolaasii; the same ASV is classified here by the RDP classifier as 
# Pseudomonas fluorescens. It is present at very low levels in ~360 samples;
# remove it
cpn.ft.filt <- cpn.ft.filt[-which(rownames(cpn.ft.filt) == "c394d6ecf2673dc653b006c53547541c"),]
cpn.tax <- cpn.tax[-which(rownames(cpn.tax) == "c394d6ecf2673dc653b006c53547541c"),]

# re-order feature table by read totals across all samples (i.e. abundance)
cpn.ft.filt <- cpn.ft.filt[order(rowSums(cpn.ft.filt),decreasing = T),]

# remove all ASVs with fewer than 500 reads across 1030 samples
cpn.ft.filt <- cpn.ft.filt[-which(rowSums(cpn.ft.filt) <500),]

# check for and remove any samples with <1,000 reads
length(which(colSums(cpn.ft.filt) <1000)) # 2 samples
cpn.ft.filt <- cpn.ft.filt[,-which(colSums(cpn.ft.filt) <1000)]

# filter taxonomy dataframe for the ASVs in the feature table
cpn.tax.filt <- cpn.tax[rownames(cpn.ft.filt),]

############################### taxonomy lookup ###############################

# extract confidence columns from df
conf <- cpn.tax.filt[,grep("conf", colnames(cpn.tax.filt))]

# for each ASV, pull the confidence score correspondong to the lowest taxonomic
# rank (i.e. how far down the taxonomic hierarchy can we classify each ASV)
conf.clas <- vector()
for(i in 1:nrow(conf)){
  conf.clas[i] <- length(which(conf[i,] >=0.8))
}

# use the above info to pull the lowest level taxonomy for each ASV
tax.low <- vector()
for(i in 1:length(conf.clas)){
  tax.low[i] <- cpn.tax.filt[i, paste0("tax",conf.clas[i])]
}

# for the full taxonomy, make a new vector containing each row taxon separated
# by two underscores and the rank initial

# pull relevant taxonomy columns
tax.all <- cpn.tax.filt[grep("tax",colnames(cpn.tax.filt))] %>% 
  dplyr::select(-c(tax1,tax2))

# add underscores and rank initials
tax.all$tax3 <- paste0("p__",tax.all$tax3)
tax.all$tax4 <- paste0("c__",tax.all$tax4)
tax.all$tax5 <- paste0("o__",tax.all$tax5)
tax.all$tax6 <- paste0("f__",tax.all$tax6)
tax.all$tax7 <- paste0("g__",tax.all$tax7)
tax.all$tax8 <- paste0("s__",tax.all$tax8)

# loop over all ASVs and pull the full taxonomy up to the lowest rank. Note that
# the 'stringr' library is needed to collapse the contents of the 'tax.all'
# dataframe subset (i.e. multiple columns of a single row) into a single string
# using 'str_flatten()' and collapsing by a semi-colon
tax.all.long <- vector()
for(i in 1:length(conf.clas)){
  j <- max(conf.clas[i])-2
  tax.all.long[i] <- tax.all[i,1:j] %>% 
    str_flatten(collapse = ";")
}

# make final df containing single taxonomic assignment, full taxonomy (single
# string) followed by each taxonomic rank (replcing NAs with 'Unclassified')

# dataframe for each taxonomic rank
final.tax <- as.data.frame(matrix(data = NA, ncol = ncol(tax.all), nrow = length(conf.clas),
                                  dimnames = list(rownames(tax.all),c("phylum","class","order","family","genus","species"))))

# add info for ranks >0.8 confidence
for(i in 1:length(conf.clas)){
  j <- max(conf.clas[i])-2
  final.tax[i,1:j] <- tax.all[i,1:j] 
}

# replace NAs with "unclassified"
for(i in 1:ncol(final.tax)){
  final.tax[,i] <- case_when(is.na(final.tax[,i]) ~ "Unclassified",
                             .default = final.tax[,i])
}

# bind single rank and all rank taxonomic assignment vectors and rearrange
final.tax <- cbind(tax=tax.low, tax.long = tax.all.long, final.tax)

# manually edit all Shigella to E.coli; not only was this done in the original
# paper (due to cpn60 seqs being incredibly similar between the two genera and
# these two taxa essentially being the same), but BLASTing the 150 bp ASVs which
# are classified to Esch/Shig vs. all cpn60 group I reference sequences from
# cpnDB shows they are ALL E.coli anyway. This explains why no ASV classified to
# Esch/Shig species have a confidence value of >0.8 (also, GTDB has collapsed 
# almost all Shigella spp. including S. flexneri, dysenteriae and sonnei into 
# E. coli anyway, so I feel justified in doing this again)

ecoli <- c("b88ede60670c601cc2bb3b58ddb4a463","01a98d33c66aa7178ad91e13973b407c",
           "5520a5205a333fc5c0d8a2936812f5e2","930e642372959f97b5aa36572afe3bcf",
           "8dac900dafb558ca6333adf04306251b","95d51701bbd86dec0acf4a28587e4d91",
           "4ca7da8e450ee976bed312c36d7d3774","98e6fb149e3825d799229e1e5aa433cf",
           "62dcda940163a1bd7d2a3cc7ffab8ee9","a542809fadabd2c2511740510acc0981",
           "025c517a5ad48252666f0c11cf8d232a","aa8641ab4e8c8f815b910ea6cd65f20f",
           "e2e03f27f92964da173274723d40263a","596308cddc3b473d84681bfc92a3754a",
           "1f77ed560f517736ed72f3ba10f398df","70224eab342c233b7ec32352209748d1",
           "075496a68bddf8e2b9d8111c03c095c5","a7a0233cbc14d2caf427d066afcbbcaa",
           "3c844156ed0c9051868390845a8926d1","828f2b2cc2e1433471efdd7e06748901",
           "af70ad87022ff8d73a771ffd59dca2b4","0e4dd0031abaf587202194568de92ec8",
           "e9f700f10039e9dd5a7821ffe4152374","b9a461a0ac5fef01cf67677dc7db14d5",
           "1a9f18ea231bf59a4d98d85bfe39583b","f2948c384c8fdfac2a5cac560339a714")

final.tax$tax.long <- case_when(rownames(final.tax) %in% ecoli ~ "p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
                                .default = final.tax$tax.long)

final.tax$tax <- case_when(rownames(final.tax) %in% ecoli ~ "Escherichia coli",
                           .default = final.tax$tax)

final.tax$genus <- case_when(rownames(final.tax) %in% ecoli ~ "g__Escherichia",
                           .default = final.tax$genus)

final.tax$species <- case_when(rownames(final.tax) %in% ecoli ~ "s__Escherichia_coli",
                           .default = final.tax$species)

final.tax$species <- gsub(" ", "_", final.tax$species)

# sanity check that ASV order is maintained in ft and tax lookup
all(rownames(cpn.ft.filt) == rownames(final.tax)) # returns TRUE

# rename ft and tax lookup and save as text
cpn60.ft <- cpn.ft.filt
cpn60.tax <- final.tax

write.table(cpn60.ft, file = "~/Documents/GitHub/usri/data/cpn60.data.txt",
            sep = "\t", quote = F, row.names = T, col.names = T)

write.table(cpn60.tax, file = "~/Documents/GitHub/usri/data/cpn60.taxa.txt",
            sep = "\t", quote = F, row.names = T, col.names = T)
