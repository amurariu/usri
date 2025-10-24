# unpermuted differential expression analysis: ALDEx2 vs. ALDEx3

# Scott Dos Santos
# Last edited: 2025-10-24

library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# set path to GH repository
repo <- "~/Documents/GitHub/usri/"

#################################### setup ####################################

# load in following from GH repo: output of ALDEx2/ALDEx3 run on unpermuted BRCA
# dataset, ensembl -> gene symbol lookup, and lists of validated/candidate 
# cancer drivers from Network of Cancer Genes v7.2 (obained 2025-10-20, source:
# http://network-cancer-genes.org/cancertypedrivers.ph )

# unpermuted ALDEx2, all gamma, brca
url.a2 <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/brca.a2.u.Rda"
tf.a2 <- tempfile(fileext = ".Rda")
download.file(url = url.a2, destfile = tf.a2, mode = "wb")
load(tf.a2)
unlink(tf.a2)

# unpermuted ALDEx3, all gamma, brca
url.a3 <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/brca.a3.u.Rda"
tf.a3 <- tempfile(fileext = ".Rda")
download.file(url = url.a3, destfile = tf.a3, mode = "wb")
load(tf.a3)
unlink(tf.a3)

# ensembl ID -> gene symbol ID lookup table
url.lookup <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/TCGA-BRCA_geneLookupAll.txt"
brca.lookup <- read.table(url.lookup, sep = "\t", header = T, quote = "", row.names = 1)

# candidate & validated genes involved in breast cancer
url.ncgbst <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/brca.NCG.breastCancer.txt"
ncg.breast <- read.table(url.ncgbst, sep = "\t", header = T, quote = "")

# candidate & validated genes involved in cancer (all types)
url.ncgall <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/brca.NCG.allCancer.txt"
ncg.all <- read.table(url.ncgall, sep = "\t", header = T, quote = "")

rm(list = c(ls(pattern = "url\\."), ls(pattern = "tf\\.")))



# # original code to generate 'brca.a2.u' and 'brca.a3.u' from scratch:
# 
# # local directory containing analysis outputs
# data <- "~/Documents/GitHub/ext_analysis/"
# 
# # lists to hold each unpermuted analysis output
# df.a2 <- list()
# df.a3 <- list()
# 
# for(i in list.files(data, pattern = "brca.data.aldex",full.names = T)){
#   
#   load(i)
# 
#   tool <- gsub("/Users/scottdossantos/Documents/GitHub/ext_analysis//brca\\.data\\.", "",
#                gsub("\\.Rda", "", i))
#   
#   if(str_split_i(ls(pattern = "brca"), pattern = "\\.", 3) == "aldex2"){
#     
#     df.a2[[tool]] <- get(ls(pattern = "brca"))$u.data
#     df.a2[[tool]]$tool <- gsub("_.*", "", tool)
#     df.a2[[tool]]$gamma <- paste0("0.", gsub("aldex.*_", "", tool))
#       
#     } else{
#       
#     df.a3[[tool]] <- get(ls(pattern = "brca"))$u.data
#     df.a3[[tool]]$tool <- gsub("_.*", "", tool)
#     df.a3[[tool]]$gamma <- paste0("0.", gsub("aldex.*_", "", tool))
#       
#     }
#     
#   rm(list = c(ls(pattern = "brca"), "tool", "i"))
#   
#   }
# 
# # collapse unpermuted results to list
# brca.a2.u <- do.call(rbind, df.a2)
# brca.a3.u <- do.call(rbind, df.a3)
# 
# # save collapsed lists to .Rda
# save(brca.a2.u, file = paste0(repo, "data/brca.a2.u.Rda"))
# save(brca.a3.u, file = paste0(repo, "data/brca.a3.u.Rda"))
# 
# # pull all Ensembl gene IDs (to extract gene symbols and descriptions from
# # BioMart at: https://useast.ensembl.org/info/data/biomart/index.html)
# # write(gsub("\\..*", "", rownames(df.a2$aldex2_0)), file = paste0(repo, "data/TCGA-BRCA_allGenes.txt"))

################################## DA: ALDEx2 ##################################

# make entity and lookup columns for each gene in the DA results
brca.a2.u$entity <- gsub("aldex2_.\\.", "", rownames(brca.a2.u))
brca.a2.u$entity <- gsub("\\..*", "", brca.a2.u$entity)

brca.a2.u$gene <- brca.lookup[match(brca.a2.u$entity, brca.lookup$Gene.stable.ID), "HGNC.symbol"]
brca.a2.u$gene.desc <- brca.lookup[match(brca.a2.u$entity, brca.lookup$Gene.stable.ID), "Gene.description"]


# at each gamma value filter aldex2 results for significant results
for(i in levels(factor(brca.a2.u$gamma))){
  assign(x = paste0("sig.a2.", i),
         value = brca.a2.u %>% 
           filter(we.eBH <0.05, gamma == i))
}

# check how many significantly different genes are NOT in the lookup
length(which(is.na(sig.a2.0.0$gene))) # 336

# identify genes which were dropped at each gamma increase
dropped.01 <- data.frame(entity = setdiff(sig.a2.0.0$entity, sig.a2.0.1$entity))
dropped.02 <- data.frame(entity = setdiff(sig.a2.0.1$entity, sig.a2.0.2$entity))
dropped.03 <- data.frame(entity = setdiff(sig.a2.0.2$entity, sig.a2.0.3$entity))
dropped.04 <- data.frame(entity = setdiff(sig.a2.0.3$entity, sig.a2.0.4$entity))
dropped.05 <- data.frame(entity = setdiff(sig.a2.0.4$entity, sig.a2.0.5$entity))

dropped.a2 <- data.frame(entity = rbind(dropped.01,dropped.02, dropped.03, dropped.04, dropped.05),
                         from = c(rep("0.1",nrow(dropped.01)), rep("0.2",nrow(dropped.02)), rep("0.3",nrow(dropped.03)), rep("0.4",nrow(dropped.04)), rep("0.5",nrow(dropped.05))))

dropped.a2$gene <- brca.lookup[match(dropped.a2$entity, brca.lookup$Gene.stable.ID), "HGNC.symbol"]
dropped.a2$desc <- brca.lookup[match(dropped.a2$entity, brca.lookup$Gene.stable.ID), "Gene.description"]

rm(list = ls(pattern = "dropped\\.0.*"))

# separate aggregated ALDEx3 results by gamma value
for(i in levels(factor(brca.a2.u$gamma))){
  gam <- gsub("0\\.", "", i)
  assign(x = paste0("a2.", gam),
         value = brca.a2.u %>% 
           filter(gamma == i))
  rm(gam,i)
}

# add ncg status (breast cancer) to ALL dropped ALDEx3 genes
dropped.a2$ncg.bst <- ncg.breast[match(dropped.a2$gene, ncg.breast$gene),2]
dropped.a2$ncg.bst <- case_when(is.na(dropped.a2$ncg.bst) ~ "n/a",
                                .default = dropped.a2$ncg.bst)

# for each dropped gene, pull dispersion, log difference for gamma = 0, and 
# qval/FDR for all gamma values. NOTE: 'std.error' and 'dispersion' do not 
# change with different gamma values, and 'estimate' *barely* changes 
dropped.a2$diff.btw <- a2.0[match(dropped.a2$entity, a2.0$entity),"diff.btw"]
dropped.a2$diff.win <- a2.0[match(dropped.a2$entity, a2.0$entity),"diff.win"]

dropped.a2$qval0.0 <- a2.0[match(dropped.a2$entity, a2.0$entity),"we.eBH"]
dropped.a2$qval0.1 <- a2.1[match(dropped.a2$entity, a2.1$entity),"we.eBH"]
dropped.a2$qval0.2 <- a2.2[match(dropped.a2$entity, a2.2$entity),"we.eBH"]
dropped.a2$qval0.3 <- a2.3[match(dropped.a2$entity, a2.3$entity),"we.eBH"]
dropped.a2$qval0.4 <- a2.4[match(dropped.a2$entity, a2.4$entity),"we.eBH"]
dropped.a2$qval0.5 <- a2.5[match(dropped.a2$entity, a2.5$entity),"we.eBH"]

# write this table to .txt, to be included in a supplementary file
# write.table(dropped.a2, file = paste0(repo, "data/SupplementaryFile1_droppedALDEx2.txt"),
#             quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

################################## DA: ALDEx3 ##################################

# make entity and lookup columns for each gene in the DA results
brca.a3.u$entity2 <- gsub("\\..*", "", brca.a2.u$entity)

brca.a3.u$gene <- brca.lookup[match(brca.a3.u$entity2, brca.lookup$Gene.stable.ID), "HGNC.symbol"]
brca.a3.u$gene.desc <- brca.lookup[match(brca.a3.u$entity2, brca.lookup$Gene.stable.ID), "Gene.description"]

# calculate dispersion (difference within groups) for ALDEx3 output: diff.win
# is calculated as std.error * sqrt(no.samples). 224 samples for BRCA
brca.a3.u <- brca.a3.u %>% 
  mutate(diff.win = std.error * sqrt(224)) %>% 
  relocate(diff.win, .after = estimate)

# at each gamma value filter aldex3 results for significant results
for(i in levels(factor(brca.a3.u$gamma))){
  assign(x = paste0("sig.a3.", i),
         value = brca.a3.u %>% 
           filter(p.val.adj <0.05, gamma == i))
}

# check how many significantly different genes are NOT in the lookup
length(which(is.na(sig.a3.0.0$gene))) # 336

# identify genes which were dropped at each gamma increase
dropped.01 <- data.frame(entity = setdiff(sig.a3.0.0$entity2, sig.a3.0.1$entity2))
dropped.02 <- data.frame(entity = setdiff(sig.a3.0.1$entity2, sig.a3.0.2$entity2))
dropped.03 <- data.frame(entity = setdiff(sig.a3.0.2$entity2, sig.a3.0.3$entity2))
dropped.04 <- data.frame(entity = setdiff(sig.a3.0.3$entity2, sig.a3.0.4$entity2))
dropped.05 <- data.frame(entity = setdiff(sig.a3.0.4$entity2, sig.a3.0.5$entity2))

dropped.a3 <- data.frame(entity = rbind(dropped.01,dropped.02, dropped.03, dropped.04, dropped.05),
                         from = c(rep("0.1",nrow(dropped.01)), rep("0.2",nrow(dropped.02)), rep("0.3",nrow(dropped.03)), rep("0.4",nrow(dropped.04)), rep("0.5",nrow(dropped.05))))

dropped.a3$gene <- brca.lookup[match(dropped.a3$entity, brca.lookup$Gene.stable.ID), "HGNC.symbol"]
dropped.a3$desc <- brca.lookup[match(dropped.a3$entity, brca.lookup$Gene.stable.ID), "Gene.description"]

rm(list = ls(pattern = "dropped\\.0.*"))

# separate aggregated ALDEx3 results by gamma value
for(i in levels(factor(brca.a3.u$gamma))){
  gam <- gsub("0\\.", "", i)
  assign(x = paste0("a3.", gam),
         value = brca.a3.u %>% 
           filter(gamma == i))
  rm(gam,i)
}

# add ncg status (breast cancer) to ALL dropped ALDEx3 genes
dropped.a3$ncg.bst <- ncg.breast[match(dropped.a3$gene, ncg.breast$gene),2]
dropped.a3$ncg.bst <- case_when(is.na(dropped.a3$ncg.bst) ~ "n/a",
                                .default = dropped.a3$ncg.bst)

# for each dropped gene, pull dispersion, log difference for gamma = 0, and 
# qval/FDR for all gamma values. NOTE: 'std.error' and 'dispersion' do not 
# change with different gamma values, and 'estimate' *barely* changes 
dropped.a3$estimate <- a3.0[match(dropped.a3$entity, a3.0$entity2),"estimate"]
dropped.a3$diff.win <- a3.0[match(dropped.a3$entity, a3.0$entity2),"diff.win"]

dropped.a3$qval0.0 <- a3.0[match(dropped.a3$entity, a3.0$entity2),"p.val.adj"]
dropped.a3$qval0.1 <- a3.1[match(dropped.a3$entity, a3.1$entity2),"p.val.adj"]
dropped.a3$qval0.2 <- a3.2[match(dropped.a3$entity, a3.2$entity2),"p.val.adj"]
dropped.a3$qval0.3 <- a3.3[match(dropped.a3$entity, a3.3$entity2),"p.val.adj"]
dropped.a3$qval0.4 <- a3.4[match(dropped.a3$entity, a3.4$entity2),"p.val.adj"]
dropped.a3$qval0.5 <- a3.5[match(dropped.a3$entity, a3.5$entity2),"p.val.adj"]

# write this table to .txt, to be included in a supplementary file
# write.table(dropped.a3, file = paste0(repo, "data/SupplementaryFile1_droppedALDEx3.txt"),
#             quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

################################ effect: ALDEx3 ################################

# calculate -log10 pval for aldex3
a3.0$qval <- -log10(a3.0$p.val.adj)

# convert Inf to reasonable value (max non-Inf is 94.5)
a3.0$qval <- case_when(a3.0$qval == Inf ~ 95, .default = a3.0$qval)

# filter gamma = 0 results by significance at each gamma value
volc.a3.sig0 <- a3.0 %>% 
  filter(p.val.adj <0.05, a3.1$p.val.adj >0.05, a3.2$p.val.adj >0.05, 
         a3.3$p.val.adj >0.05, a3.4$p.val.adj >0.05, a3.5$p.val.adj >0.05)

volc.a3.sig1 <- a3.0 %>% 
  filter(a3.1$p.val.adj <0.05, a3.2$p.val.adj >0.05, a3.3$p.val.adj >0.05,
         a3.4$p.val.adj >0.05, a3.5$p.val.adj >0.05)

volc.a3.sig2 <- a3.0 %>% 
  filter(a3.2$p.val.adj <0.05, a3.3$p.val.adj >0.05, a3.4$p.val.adj >0.05, 
         a3.5$p.val.adj >0.05)

volc.a3.sig3 <- a3.0 %>% 
  filter(a3.3$p.val.adj <0.05, a3.4$p.val.adj >0.05, a3.5$p.val.adj >0.05)

volc.a3.sig4 <- a3.0 %>% 
  filter(a3.4$p.val.adj <0.05, a3.5$p.val.adj >0.05)

volc.a3.sig5 <- a3.0 %>% 
  filter(a3.5$p.val.adj <0.05)

volc.a3.ns <- a3.0 %>% 
  filter(p.val.adj >0.05)

# set colours for gamma values for ALDEx3 as in FDR/TPR figure
cols.volcano <- c("#E6E6FA","#CBE3FC","#B0E2FF","#66B3F6","#2171B5","#08306B",
                  brewer.pal(n = 9, name = "YlOrRd")[c(1,3,5,6,8,9)])

leg.cols.a3 <- c("0" = cols.volcano[7], "0.1" = cols.volcano[8], "0.2" = cols.volcano[9],
                 "0.3" = cols.volcano[10], "0.4" = cols.volcano[11], "0.5" = cols.volcano[12])

# set strip title
volc.a3.ns$title <- "ALDEx3: unpermuted BRCA dataset"

# png(paste0(repo, "figures/fig_unpermutedDA_ALDEx3.png"),
#     units = "in", height = 6, width = 10, res = 600)

volc.a3 <- ggplot(data = volc.a3.ns, aes(x = estimate, y = qval))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = volc.a3.sig0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+  
  geom_point(data = volc.a3.sig5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+  # FDR threshold
  scale_y_continuous(limits = c(0,97.5), expand = c(0.005,0.75))+
  scale_x_continuous(limits = c(-7.55,8.25), expand = c(0.005, 0.005))+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a3)+
  labs(x = "Log difference between groups", y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10))

volc.a3

# dev.off()

# edited version with no y axis

volc.a3.edit <- ggplot(data = volc.a3.ns, aes(x = estimate, y = qval))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = volc.a3.sig0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a3.sig4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+  
  geom_point(data = volc.a3.sig5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+  # FDR threshold
  scale_y_continuous(limits = c(0,97.5), expand = c(0.005,0.75))+
  scale_x_continuous(limits = c(-7.55,8.25), expand = c(0.005, 0.005))+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a3)+
  labs(x = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10),
        axis.title.y = element_blank())

volc.a3.edit

################################ effect: ALDEx2 ################################

# calculate -log10 pval for aldex3
a2.0$qval <- -log10(a2.0$we.eBH)

# convert Inf to reasonable value (max non-Inf is 91.3)
a2.0$qval <- case_when(a2.0$qval == Inf ~ 92.5, .default = a2.0$qval)

# filter gamma = 0 results by significance at each gamma value
volc.a2.sig0 <- a2.0 %>% 
  filter(we.eBH <0.05, a2.1$we.eBH >0.05, a2.2$we.eBH >0.05, 
         a2.3$we.eBH >0.05, a2.4$we.eBH >0.05, a2.5$we.eBH >0.05)

volc.a2.sig1 <- a2.0 %>% 
  filter(a2.1$we.eBH <0.05, a2.2$we.eBH >0.05, a2.3$we.eBH >0.05,
         a2.4$we.eBH >0.05, a2.5$we.eBH >0.05)

volc.a2.sig2 <- a2.0 %>% 
  filter(a2.2$we.eBH <0.05, a2.3$we.eBH >0.05, a2.4$we.eBH >0.05, 
         a2.5$we.eBH >0.05)

volc.a2.sig3 <- a2.0 %>% 
  filter(a2.3$we.eBH <0.05, a2.4$we.eBH >0.05, a2.5$we.eBH >0.05)

volc.a2.sig4 <- a2.0 %>% 
  filter(a2.4$we.eBH <0.05, a2.5$we.eBH >0.05)

volc.a2.sig5 <- a2.0 %>% 
  filter(a2.5$we.eBH <0.05)

volc.a2.ns <- a2.0 %>% 
  filter(we.eBH >0.05)

leg.cols.a2 <- c("0" = cols.volcano[1], "0.1" = cols.volcano[2], "0.2" = cols.volcano[3],
                 "0.3" = cols.volcano[4], "0.4" = cols.volcano[5], "0.5" = cols.volcano[6])

# set strip title
volc.a2.ns$title <- "ALDEx2: unpermuted BRCA dataset"

# png(paste0(repo, "figures/fig_unpermutedDA_ALDEx2.png"),
#     units = "in", height = 6, width = 10, res = 600)

volc.a2 <- ggplot(data = volc.a2.ns, aes(x = diff.btw, y = qval))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = volc.a2.sig0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a2.sig1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a2.sig2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a2.sig3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a2.sig4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = volc.a2.sig5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+  # FDR threshold
  scale_y_continuous(limits = c(0,97.5), expand = c(0.005,0.75))+
  scale_x_continuous(limits = c(-7.55,8.35), expand = c(0.005, 0.005))+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a2)+
  labs(x = "Log difference between groups", y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10))

volc.a2

# dev.off()

# plot both A2 and A3 side by side
# png(paste0(repo, "figures/fig_unpermutedDA_ALDExBoth.png"),
#     units = "in", height = 6, width = 12, res = 600)

volc.a2 | volc.a3.edit

# dev.off()
