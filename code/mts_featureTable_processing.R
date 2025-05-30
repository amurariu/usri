# feature table processing: vaginal metatranscriptome

# Scott Dos Santos
# 2025-05-30

# takes the raw feature table and relevant lookup tables for a vaginal
# metatranscriptome dataset (see: https://doi.org/10.1186/s40168-024-01992-w) 
# and produces a ready-to-analyse feature table of counts per KEGG orthology 
# term per sample (i.e. aggregated at the functional level)

#################################### setup ####################################

# load required packages
library(CoDaSeq)
library(dplyr)

# locations of metatranscriptome feature table (output from VIRGO), lookup table
# of 1) vNumbers to kegg orthology term; 2) vNumbers to kegg pathway; 3) KO term
# to pathway (non-redundant, manually edited)
loc.ft <- "https://raw.githubusercontent.com/amurariu/usri/main/data/mts_virginiaRawFeatureTable.txt.zip"
loc.vnum.to.ko <- "https://raw.githubusercontent.com/amurariu/usri/main/data/mts_virgo_vnumKO.txt"
loc.vnum.to.path <- "https://raw.githubusercontent.com/amurariu/usri/main/data/mts_virgo_vnumPath.txt"
loc.ko.to.path <- "https://raw.githubusercontent.com/amurariu/usri/main/data/mts_kegg_KOpath.txt"

# create temporary file for zipped feature table from GitHub
temp.zip <- tempfile(fileext = ".zip")

# download from GH into temp file
download.file(loc.ft, destfile = temp.zip, mode = "wb")

# read in above dataframes 
df <- read.table(unz(description = temp.zip, filename = "t2t_summary.virginia.NR.abundance.txt"),
                 header = T, row.names = 1, quote = '', sep = "\t")

ko <- read.table(file = loc.vnum.to.ko, header = F, sep = "\t", quote = "", row.names = 1)
path <- read.table(file = loc.vnum.to.path, header = F, sep = "\t", quote = "", row.names = 1)
kegg.nr <- read.table(file = loc.ko.to.path, header = T, sep = "\t", quote = "", row.names = 1)

############################## filtering lookups ##############################

# remove '_1' from feature table column names
colnames(df) <- gsub("_1", "", colnames(df))

# filter VIRGO KO lookup table by vNumbers in the dataset (takes about 10 mins)
ko.df <- ko[rownames(df),1, drop = F]

# remove NAs from vnumbers that don't have an assigned KO (approx. half)
ko.df <- ko.df %>%
  filter(!is.na(V2))

# filter VIRGO pathway lookup table by vNumbers in the dataset (takes ages!)
path.df <- path[rownames(df),]

# remove NAs from vnumbers that don't have an assigned pathway (approx. 20 mins)
path.df <- path.df %>%
  filter(!is.na(V2))

# add column to merge on KO
ko.df$id <- rownames(ko.df)
path.df$id <- rownames(path.df)

# merge tables on ID
func.df <- ko.df %>% 
  left_join(y = path.df, by = "id")

# edit column names, then set rownames to KOs
colnames(func.df) <- c ("ko","id","map","pathway")
rownames(func.df) <- func.df$id
func.df$id <- NULL

# edit NAs for cases where genes had a KO term but no assigned pathway
func.df$map <- case_when(is.na(func.df$map) ~ "NotAssigned", .default = func.df$map)
func.df$pathway <- case_when(is.na(func.df$pathway) ~ "NotAssigned", .default = func.df$pathway)

################################## KO editing ##################################

# check if all KOs present in mts data are present in lookup
all(levels(factor(func.df$ko)) %in% rownames(kegg.nr))

# check which ones are NOT present in master lookup
ko.notPresent <- levels(factor(func.df$ko))[which(!levels(factor(func.df$ko)) %in% rownames(kegg.nr))]

# quick google search of the first one suggests a KEGG name change issue:
# K00050 is not found on kegg, but the page redirects to 'K12972'. Another page
# shows K00050 having the same function as K12972, suggesting k00050 has been
# deprecated for some reason. This holds for other KOs so  manually edit missing
# KOs in func.df
func.df$ko <- case_when(func.df$ko == "K00050" ~ "K12972",
                        func.df$ko == "K00210" ~ "K04517",
                        func.df$ko == "K00329" ~ "K00330",
                        func.df$ko == "K00356" ~ "K03885",
                        func.df$ko == "K00536" ~ "K03839",
                        func.df$ko == "K01234" ~ "K01208",
                        func.df$ko == "K01263" ~ "K11140",
                        func.df$ko == "K01541" ~ "K01542",
                        func.df$ko == "K01572" ~ "K01571",
                        func.df$ko == "K01786" ~ "K03077",
                        func.df$ko == "K02428" ~ "K01519",
                        func.df$ko == "K02802" ~ "K02803",
                        func.df$ko == "K02808" ~ "K02810",
                        func.df$ko == "K03080" ~ "K03077",
                        func.df$ko == "K03081" ~ "K03078",
                        func.df$ko == "K03082" ~ "K03079",
                        func.df$ko == "K03815" ~ "K03783",
                        .default = func.df$ko)

# pull pathway and map values for each VIRGO gene & KO term from the manually 
# curated, non-redundant KO -> pathway lookup, then correct for deprecated KOs
func.df$pathway <- kegg.nr[func.df$ko, "path.desc"]
func.df$map <- kegg.nr[func.df$ko, "path.id"]
func.df$pathway <- case_when(is.na(func.df$pathway) ~ "NotAssigned", .default = func.df$pathway)

# pull ko descriptions for each VIRGO gene & KO term from the non-redundant,
# manually curated KO -> pathway lookup, then correct for deprecated KOs
func.df$product <- kegg.nr[func.df$ko, "ko.desc"]
func.df$product <- case_when(is.na(func.df$product) ~ "NotAssigned", .default = func.df$product)

# rearrange columns for convenience
func.df <- func.df %>% 
  relocate(map, .after = product)

# editing of specific KO terms to correct pathway as per metatranscriptome
# paper (all terms confirmed by manual search on KEGG)
func.df$pathway <-case_when(func.df$ko == "K01580" ~ "ABC transporters",
                            func.df$ko == "K01580" ~ "Acid survival",
                            func.df$ko %in% c("K00812","K11358","K01915","K01775","K01776","K00811") ~ "Alanine aspartate and glutamate metabolism",
                            func.df$ko == "K14195" ~ "Bacterial adherence",
                            func.df$ko == "K06595" ~ "Bacterial chemotaxis",
                            func.df$ko == "K11936" ~ "Exopolysaccharide biosynthesis",
                            func.df$ko == "K14534" ~ "Butanoate metabolism",
                            func.df$ko %in% c("K03367","K03739","K03740","K14188","K14205") ~ "Cationic antimicrobial peptide (CAMP) resistance",
                            func.df$ko %in% c("K03531","K03588","K03589","K03590","K13581") ~ "Cell division",
                            func.df$ko %in% c("K00024","K00174","K00175","K00239","K00240","K00241","K01676","K01681","K01847","K01848","K01849","K01902","K01903","K01966","K03737","K05606","K13788","K00031","K00177","K00625","K00925","K01595","K01895","K05884","K01679") ~ "Citrate cycle (TCA cycle)",
                            func.df$ko %in% c("K01358","K03544") ~ "Clp-dependent proteolysis",
                            func.df$ko == "K11031" ~ "Cytolysis",
                            func.df$ko %in% c("K02313","K02314") ~ "DNA replication",
                            func.df$ko == "K00864" ~ "Glycerolipid metabolism",
                            func.df$ko == "K00981" ~ "Glycerophospholipid metabolism",
                            func.df$ko %in% c("K02446","K04041","K00134","K01803","K00873","K00927","K01610","K01624","K00850","K01834","") ~ "Glycolysis / Gluconeogenesis",
                            func.df$ko == "K00018" ~ "Glyoxylate and dicarboxylate metabolism",
                            func.df$ko %in% c("K00058","K00600","K01079") ~ "Glycine serine and threonine metabolism",
                            func.df$ko == "K01715" ~ "Fatty acid biosynthesis",
                            func.df$ko %in% c("K02405", "K02406", "K03092","K06603","K13626") ~ "Flagellar assembly",
                            func.df$ko == "K00817" ~ "Histidine metabolism",
                            func.df$ko %in% c("K11749","K02012","K02217","K07243") ~ "Iron acquisition",
                            func.df$ko == "K01338" ~ "Lon protease",
                            func.df$ko == "K01582" ~ "Lys metabolism",
                            func.df$ko %in% c("K05299","K00123","K00124","K00127") ~ "Methanoate metabolism",
                            func.df$ko %in% c("K04079","K04043","K04077") ~ "Molecular chaperones",
                            func.df$ko %in% c("K00297","K01491","K01938") ~ "One carbon pool by folate",
                            func.df$ko %in% c("K02117","K02118","K02119","K02120","K02121","K02122","K02123","K02124") ~ "Oxidative phosphorylation",
                            func.df$ko == "K00832" ~ "Phenylalanine tyrosine and tryptophan metabolism",
                            func.df$ko == "K01195" ~ "Pentose and glucuronate interconversions",
                            func.df$ko %in% c("K01621","K01632","K01783","K01807","K01808","K00615") ~ "Pentose phosphate pathway",
                            func.df$ko == "K02563" ~ "Peptidoglycan biosynthesis",
                            func.df$ko == "K01488" ~ "Purine metabolism",
                            func.df$ko == "K00757" ~ "Pyrimidine metabolism",
                            func.df$ko %in% c("K00029","K01006","K01007") ~ "Pyruvate metabolism",
                            func.df$ko == "K01661" ~ "Quinone biosynthesis",
                            func.df$ko %in% c("K14051","K07173") ~ "Quorum sensing",
                            func.df$ko == "K01186" ~ "Sialidase activity",
                            func.df$ko %in% c("K03313","K03315") ~ "Sodium-proton antiport",
                            func.df$ko == "K00688" ~ "Starch and sucrose metabolism",
                            func.df$ko == "K04564" ~ "ROS degradation",
                            func.df$ko == "K00869" ~ "Terpenoid biosynthesis",
                            func.df$ko == "K02358" ~ "Translation",
                            func.df$ko == "K02040" ~ "Two-component system",
                            func.df$ko == "K08303" ~ "Type I collagen degradation",
                            func.df$ko == "K00532" ~ "Unknown",
                            func.df$ko == "K01428" ~ "Urease activity",
                            func.df$ko %in% c("K03183","K03688") ~ "Ubiquinone biosynthesis",
                            func.df$ko %in% c("K02548","K02548") ~ "Vitamin K2 biosynthesis",
                            .default = func.df$pathway)

# even further editing of specific KO terms to correct pathway (discovered
# after analysing data downstream- these were all 'Unknown'; references to
# support changes are in the code of the GH repository for the metatranscriptome
# paper)
func.df$pathway <- case_when(func.df$ko %in% c("K00435","K01772") ~ "Porphyrin metabolism",
                             func.df$ko == "K07777" ~ "Two-component system",
                             func.df$ko == "K01005" ~ "Teichoic acid biosynthesis",
                             func.df$ko == "K00903" ~ "Exopolysaccharide biosynthesis",
                             func.df$ko %in% c("K01356","K03631","K03700") ~ "DNA damage repair",
                             func.df$ko == "K02055" ~ "Putrescine/Spermidine transport",
                             func.df$ko %in% c("K02068","K02069") ~ "Iron homeostasis",
                             func.df$ko %in% c("K02237","K02244","K02248") ~ "Extracellular DNA uptake",
                             func.df$ko %in% c("K02444","K03436") ~ "Fructose metabolism",
                             func.df$ko %in% c("K03282","K03442") ~ "Osmotic stress adaptation",
                             func.df$ko == "K02530" ~ "Lactose metabolism",
                             func.df$ko == "K02859" ~ "Riboflavin metabolism",
                             func.df$ko == "K03284" ~ "Magnesium transport",
                             func.df$ko %in% c("K03297","K08161") ~ "Multidrug resistance - efflux pump",
                             func.df$ko == "K03303" ~ "Lactate metabolism",
                             func.df$ko == "K03322" ~ "Manganese acquisition",
                             func.df$ko %in% c("K03484", "K05992") ~ "Starch and sucrose metabolism",
                             func.df$ko == "K04075" ~ "Aminoacyl-tRNA biosynthesis",
                             func.df$ko == "K07043" ~ "Pyrimidine metabolism",
                             func.df$ko == "K01595" ~ "Pyruvate metabolism",
                             func.df$ko == "K01066" ~ "Unknown",
                             .default = func.df$pathway)

# generate ko vector and collapse read dataframe by KO term
ko.vec <- func.df[rownames(df),"ko"]
names(ko.vec) <- rownames(df)

mts.data <- aggregate(df, by = list(ko.vec), FUN = sum)
rownames(mts.data) <- mts.data$Group.1
mts.data$Group.1 <- NULL


# remove eukaryotic KOs from table 
mts.data <- mts.data[-which(grepl(paste("K03260","K06237","K12373","K00863",
                                        "K13993","K03648","K01408","K00599", sep = "|"),
                                  rownames(mts.data))),]

# filter gene-level dataset and functional lookup table using original parameters
df.filt <- codaSeq.filter(df, samples.by.row = FALSE,
                          min.occurrence=0.30, min.prop=0.00005)

func.df.filt <- func.df[rownames(df.filt),]
rownames(func.df.filt) <- rownames(df.filt)

# edit NAs to 'NotAssigned' for genes lacking functional information in VIRGO
func.df.filt$ko <-case_when(is.na(func.df.filt$ko) ~ "NotAssigned", .default = func.df.filt$ko)
func.df.filt$pathway <-case_when(is.na(func.df.filt$pathway) ~ "NotAssigned", .default = func.df.filt$pathway)
func.df.filt$product <-case_when(is.na(func.df.filt$product) ~ "NotAssigned", .default = func.df.filt$product)
func.df.filt$map <-case_when(is.na(func.df.filt$map) ~ "NotAssigned", .default = func.df.filt$map)

# generate ko vector and collapse filtered read dataframe by KO term
ko.vec.filt <- func.df.filt[rownames(df.filt), "ko"]
names(ko.vec.filt) <- rownames(df.filt)

mts.data.filt <- aggregate(df.filt, by = list(ko.vec.filt), FUN = sum)
rownames(mts.data.filt) <- mts.data.filt$Group.1
mts.data.filt$Group.1 <- NULL

# remove eukaryotic KOs (and 'NotAssigned') from table 
mts.data.filt <- mts.data.filt[-which(grepl(paste("K03260","K06237","K12373","K00863",
                                                  "K13993","K03648","K01408","K00599",
                                                  "NotAssigned", sep = "|"),
                                            rownames(mts.data.filt))),]

# save feature tables and vNumber -> KO -> pathway lookup as tab delimited text
write.table(x = mts.data, file = "~/Documents/GitHub/usri/data/mts.data.txt",
            quote = F, sep = "\t", row.names = T, col.names = T)

write.table(x = func.df, file = "~/Documents/GitHub/usri/data/mts.func.txt",
            quote = F, sep = "\t", row.names = T, col.names = T)

write.table(x = mts.data.filt, file = "~/Documents/GitHub/usri/data/mts.data.filt.txt",
            quote = F, sep = "\t", row.names = T, col.names = T)

write.table(x = func.df.filt, file = "~/Documents/GitHub/usri/data/mts.func.filt.txt",
            quote = F, sep = "\t", row.names = T, col.names = T)
