# processing tiyani data for differential abundance scale benchmarking

# Scott Dos Santos
# Last edited: 8th May 2026

# want to extract all pairwise comparisons from tiyani data (16S comparison of
# gut microbiome from wide age range of people)

library(dplyr)

#################################### setup ####################################

# get urls for github data
url.count <- "https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/tiyani_counts_orig.txt"
url.meta <- "https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/tiyani_meta_minimal.txt"
url.tax <- "https://raw.githubusercontent.com/amurariu/usri/refs/heads/main/data/tiyani_tax_orig.txt"

# load count table, metadata and taxonomy lookup (taken from Greg's github
# at: https://github.com/ggloor/datasets/clean_[data | meta | tax].txt)

# NOTE: several sample IDs start with numbers, which would be prefixed with 'X'
#       if 'check.names' is not set to FALSE
tiyani.ft <- read.table(url.count, header = TRUE, sep = "\t",
                        quote = "", row.names = 1, check.names = F)

tiyani.meta <- read.table(url.meta, header = TRUE, sep = "\t",
                          quote = "", row.names = 1)

tiyani.tax <- read.table(url.tax, header = TRUE, sep = "\t",
                         quote = "", row.names = 1)

# replace row names
rownames(tiyani.ft) <- rownames(tiyani.tax) <- paste0("OTU_", 1:nrow(tiyani.tax))

# check for equality of metadata row names and count table column names
all(colnames(tiyani.ft) == rownames(tiyani.meta)) # returns TRUE

# remove young soliders
tiyani.meta <- tiyani.meta %>% 
  filter(Group != "young soldier")

# add age range column
tiyani.meta$age.range <- case_when(tiyani.meta$Group == "Centenarians" ~ ">94",
                                   tiyani.meta$Group == "elder" ~ "60-79",
                                   tiyani.meta$Group == "kindergarten" ~ "3-6",
                                   tiyani.meta$Group == "mid_age" ~ "30-50",
                                   tiyani.meta$Group == "mid_school" ~ "13-14",
                                   tiyani.meta$Group == "youth" ~ "19-24",
                                   tiyani.meta$Group == "Pupils" ~ "8-12")

# edit groups to look nicer
levels(factor(tiyani.meta$Group))
tiyani.meta$Group <- case_when(tiyani.meta$Group == "Centenarians" ~ "Centenarians",
                               tiyani.meta$Group == "elder" ~ "Elderly",
                               tiyani.meta$Group == "kindergarten" ~ "Kindergarten",
                               tiyani.meta$Group == "mid_age" ~ "Middle_aged",
                               tiyani.meta$Group == "mid_school" ~ "Middle_school",
                               tiyani.meta$Group == "youth" ~ "Young_adults",
                               tiyani.meta$Group == "Pupils" ~ "Primary_school")

# group ages:
# _____________________________
#     GROUP               AGES
# _____________________________
#     Kindergarten         3-6
#     Primary_school      8-12
#     Middle_school      13-14
#     Young_adults       19-24
#     Middle_aged        30-50
#     Elderly            60-79
#     Centenarians         >94
# _____________________________

# set metadata column names to lowercase
colnames(tiyani.meta) <- c("group", "age", "health", "age.range")

#################################### combos ####################################

# get all unique pairs of groups for comparisons
pairs <- levels(factor(tiyani.meta$group))
combos <- combn(x = pairs, m = 2)

# make dataframe of groups and abbreviated names for assigning object IDs
id <- c("cent", "eld", "kind", "mage", "msch", "psch", "young")
id.df <- data.frame(name = pairs, short = id)

# set directory to write .Rda files to
outdir <- "~/Documents/GitHub/usri/data/tiyani_pairs/"

# loop over original dataset and filter for groups in the unique pairs to
# produce metadata, taxonomy and counts for all 21 pairwise combinations
if(length(list.files(outdir, pattern = "tiyani")) != 63){
  for(i in 1:ncol(combos)){
    
    pair <- combos[,i]
    pair.id <- c(id.df$short[match(pair[1], id.df$name)],
                 id.df$short[match(pair[2], id.df$name)])
    
    meta <- tiyani.meta[which(tiyani.meta$group %in% pair),]
    
    read <- tiyani.ft[,rownames(meta)]
    read <- read[-which(rowSums(read) == 0),]
    
    tax <- tiyani.tax[rownames(read),]
    
    name.ft <- paste0(pair.id[1], ".", pair.id[2], ".ft")
    name.meta <- paste0(pair.id[1], ".", pair.id[2], ".meta")
    name.tax <- paste0(pair.id[1], ".", pair.id[2], ".tax")
    
    path.ft <- paste0(wd, "tiyani_", paste0(pair.id[1], ".", pair.id[2], ".ft"), ".Rda")
    path.meta <- paste0(wd, "tiyani_", paste0(pair.id[1], ".", pair.id[2], ".meta"), ".Rda")
    path.tax <- paste0(wd, "tiyani_", paste0(pair.id[1], ".", pair.id[2], ".tax"), ".Rda")
    
    assign(x = name.ft, value = get("read"))
    assign(x = name.meta, value = get("meta"))
    assign(x = name.tax, value = get("tax"))
    
    save(list = name.ft, file = path.ft)
    save(list = name.meta, file = path.meta)
    save(list = name.tax, file = path.tax)
    
  }
}
