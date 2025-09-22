# density plots and strip chart for example dataset: ALDEx2 vs. ALDEx3

# Andreea Murariu & Scott Dos Santos
# Last edited: 2025-08-19

library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# set paths to directory containing all analysis objects/repo
data <- "~/Documents/GitHub/ext_analysis/"
repo <- "~/Documents/GitHub/usri/"

############################ summarise min diff.btw ############################

# source function for extracting relevant data
source('code/summary.data.fun.R')

# load in analysis data and get summary objects for all objects (or load them 
# if they already exist as .Rda files) 

# NOTE: Same deal for this section as in figure 1 sensitivity/fdr plots: I had 
#       to run the for loops below individually to create the 'ss.' objects 
#       for each dataset to avoid memory issues on local machine. Code is set
#       up to load all 'ss.' objects in if they all exist in the corresponding
#       directory.

if(length(list.files(paste0(repo,"analysis/summarystats"), pattern = "ss")) != 10){
  
  # load in analysis results (takes a few minutes!)
  for(i in list.files(data, pattern = "sccyto.*aldex")){
    load(paste0(data, i))
  }
  
  # extract dataset and tool names from input file
  i <- ls(pattern = "sccyto")[1]
  
  dataset <- str_split(i, "\\." , 3)[[1]][1]
  gamma <- "data_"
  tool <- "aldex"
    
  # run 'summary statistics' on analysis objects to extract the (absolute)
  # minimum difference between group values that is needed for a significant
  # BH-corrected P value at different gamma values (absolute as diff.btw
  # will be positive or negative depending on the direction of the difference).
  # max = positive values, min = negative values, abs = absolute values of
  # min column followed by absolute values of max column.
  assign(x = paste0("ss.", dataset),
         sum.fun(aldex2_0 = get(paste0(dataset, ".", gamma, 0, ".", tool, 2)),
                 aldex2_1 = get(paste0(dataset, ".", gamma, 1, ".", tool, 2)),
                 aldex2_2 = get(paste0(dataset, ".", gamma, 2, ".", tool, 2)),
                 aldex2_3 = get(paste0(dataset, ".", gamma, 3, ".", tool, 2)),
                 aldex2_4 = get(paste0(dataset, ".", gamma, 4, ".", tool, 2)),
                 aldex2_5 = get(paste0(dataset, ".", gamma, 5, ".", tool, 2)),
                 aldex3_0 = get(paste0(dataset, ".", gamma, 0, ".", tool, 3)),
                 aldex3_1 = get(paste0(dataset, ".", gamma, 1, ".", tool, 3)),
                 aldex3_2 = get(paste0(dataset, ".", gamma, 2, ".", tool, 3)),
                 aldex3_3 = get(paste0(dataset, ".", gamma, 3, ".", tool, 3)),
                 aldex3_4 = get(paste0(dataset, ".", gamma, 4, ".", tool, 3)),
                 aldex3_5 = get(paste0(dataset, ".", gamma, 5, ".", tool, 3))))
  
  
  # save summary list objects as .Rda
  # NOTE: for MTS dataset, iteration 73 of ALDEx2 at gamma = 0.5 had NO
  # significant P values for the negative diff.btw features and so returns
  # -Inf for the 73rd value in 'max' column and 173rd value in 'abs' column.
  # Same issue present in single cell dataset: 28 Inf values (3 in ALDEx2 - 0.3,
  # 7 in ALDEx2 - 0.4, 13 in ALDEx2 - 0.5, 1 in ALDEx3 - 0.4, 4 in ALDEx3 - 0.5)
    save(list = paste0("ss.", dataset),
         file = paste0(repo,"analysis/summarystats/ss.", dataset, ".Rda"))
  
  # finally, remove analysis objects from environment
  for(i in ls(pattern = "^sccyto")){
    rm(list = i)
  }
    
} else{
  
  # load in summary objects from .Rda
  for(i in list.files(paste0(repo,"analysis/summarystats"), pattern = "^ss")){
    load(paste0(repo, "analysis/summarystats/", i))
  }
}

############################# data transformation #############################

tmplist <- list()

for(i in ls(pattern = "^ss.")){
  # pull list
  df <- get(i)
  
  # make a temp dataframe with 2400 rows (200 x 12 tool/gamma combos) and 3 cols:
  # abs values, dataset and tool
  tempdf <- data.frame(matrix(data = NA, nrow = 2400, ncol = 3))
  colnames(tempdf) <- c("abs", "dataset", "tool")
  
  start <- 1
  end <- 200
  
  # loop over list elements for each dataset and pull out abs values from each
  # tool/gamma combo; stick in temp df & increment start/end vals
  for(j in 1:length(df)){
    
    tempdf[(start:end), 1] <- df[[j]][,3]
    name <- gsub("minmax", "aldex", names(df)[j])
    tempdf[start:end, 3] <- name
    
    start <- start + 200
    end <- end + 200
  }
  
  # fill in current dataset
  tempdf[,2] <- gsub("^ss\\.", "", i)
  
  # add temp dataframe of abs/dataset/tool values to list
  tmplist[[i]] <- tempdf
  
  # remove temp objects
  rm(df, end, j, name, start, tempdf)
}

# collapse list to df, then remove rownames, convert 'coeff' to numeric and
# change order of columns
plot.df <- do.call(rbind, tmplist)

# remove problematic Inf rows for MTS/sccyto datasets
length(which(plot.df$abs == Inf)) # 1 for mts, 28 for sccyto
plot.df <- plot.df[-grep(Inf, plot.df$abs),]

# split tool into tool and gamma (character and numeric)
plot.df$gamma.char <- as.character(gsub("aldex._", "0.", plot.df$tool))
plot.df$gamma.char <- as.character(gsub("0\\.0", "0", plot.df$gamma.char))
plot.df$gamma.num <- as.numeric(gsub("aldex._", "0.", plot.df$tool))
plot.df$tool <- gsub("_.$","", plot.df$tool)

# convert tools to correct case
plot.df$tool <- gsub("aldex","ALDEx", plot.df$tool)

# add 12 colours for scale values: 6 for ALDEx2, 6 for ALDEx3
#(note for future Scott: you selected these colours by running the following:
# colorRampPalette(colors = c("lavender", "lightskyblue1", "dodgerblue2", "#2171B5", "#08306B"), space = "rgb")(9))
# brewer.pal(n = 9, name = "YlOrRd")
# and choosing the ones with the most difference between them
cols.density <- c("#E6E6FA","#CBE3FC","#B0E2FF","#66B3F6","#2171B5","#08306B",
                  brewer.pal(n = 9, name = "YlOrRd")[c(1,3,5,6,8,9)])

#################### Density Plot Panel ######################

# filter for immuno dataset
plot.filter <- plot.df %>% 
  filter(dataset == "immuno")

# add colour values to filtered dataset
plot.filter$denscol <- case_when(plot.filter$gamma.char == "0" & plot.filter$tool == "ALDEx2" ~ cols.density[1],
                                 plot.filter$gamma.char == "0.1" & plot.filter$tool == "ALDEx2" ~ cols.density[2],
                                 plot.filter$gamma.char == "0.2" & plot.filter$tool == "ALDEx2" ~ cols.density[3],
                                 plot.filter$gamma.char == "0.3" & plot.filter$tool == "ALDEx2" ~ cols.density[4],
                                 plot.filter$gamma.char == "0.4" & plot.filter$tool == "ALDEx2" ~ cols.density[5],
                                 plot.filter$gamma.char == "0.5" & plot.filter$tool == "ALDEx2" ~ cols.density[6],
                                 plot.filter$gamma.char == "0" & plot.filter$tool == "ALDEx3" ~ cols.density[7],
                                 plot.filter$gamma.char == "0.1" & plot.filter$tool == "ALDEx3" ~ cols.density[8],
                                 plot.filter$gamma.char == "0.2" & plot.filter$tool == "ALDEx3" ~ cols.density[9],
                                 plot.filter$gamma.char == "0.3" & plot.filter$tool == "ALDEx3" ~ cols.density[10],
                                 plot.filter$gamma.char == "0.4" & plot.filter$tool == "ALDEx3" ~ cols.density[11],
                                 plot.filter$gamma.char == "0.5" & plot.filter$tool == "ALDEx3" ~ cols.density[12])

# convert colours to factor based on tool & gamma order
plot.filter$denscol <- factor(plot.filter$denscol, levels = cols.density)

# make vector of scale values for legend
leg.vals <- c(paste0("A2: ", c(0,0.1,0.2,0.3,0.4,0.5)),
              paste0("A3: ", c(0,0.1,0.2,0.3,0.4,0.5)))


# get group means for given dataset and add colours
group_means.filter<- plot.filter %>%
  group_by(tool, gamma.char) %>%
  summarise(mean_value = mean(abs))

group_means.filter$linecols <- cols.density

# plot example density plot for given dataset
# png(paste0(repo,"figures/fig3_densityPlotsALDEx.png"),
#     units = "in", height = 4, width = 6, res = 600)

p1 <- plot.filter %>% 
  ggplot(aes(x = abs)) + 
  geom_density(aes(fill = denscol), linewidth = 0.3, alpha=0.7)+
  scale_fill_identity(name = "Scale (\u03b3)", guide = "legend", labels = leg.vals)+
  geom_vline(data = group_means.filter, aes(xintercept = mean_value, colour=linecols), linetype="dashed", linewidth=0.3) +
  scale_colour_identity()+
  xlab('Minimum absolute log difference for significance')+ ylab('Density')+
  scale_x_continuous(limits = c(0,1.99), expand = c(0,0))+
  scale_y_continuous(limits = c(0,11), expand = c(0,0))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  theme_bw()+
  facet_wrap(~tool)+
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10),
        strip.text = element_text(face = "bold"), legend.box.spacing = unit(0.15, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        legend.key.spacing.y = unit(0.1, "cm"), legend.key.spacing.x = unit(0.2, "cm"),
        legend.key.size = unit(0.3, "cm"), legend.margin = margin(0,0,0,0.0,"cm"),
        legend.position = "top")

p1

# dev.off()

# calculate group means for ALL datasets, then filter for all human-derived, 
# bulk transcriptome datasets, and other datasets, separately
group_means.all <- plot.df %>%
  group_by(tool, dataset, gamma.num) %>%
  summarise(mean_value = mean(abs))

# also add column ensuring these data points are represented in legend
group_means.hum <- group_means.all %>% 
  filter(dataset %in% c("brca","immuno","kirc","lihc","luad","prad","thca")) %>% 
  mutate(shp = "TCGA transcriptomes")

# also add column for strip title and replace dataset abbreviations with full
# labels
group_means.oth <- group_means.all %>% 
  filter(!dataset %in% c("brca","immuno","kirc","lihc","luad","prad","thca")) %>% 
  mutate(title = "Minimum absolute difference between groups for P <0.05")

group_means.oth$dataset <- case_when(group_means.oth$dataset == "mts" ~ "Vaginal metatranscriptome",
                                     group_means.oth$dataset == "sccyto" ~ "Single-cell transcriptome",
                                     group_means.oth$dataset == "yeast" ~ "Yeast transcriptome")

group_means.oth2 <- plot.df %>% 
  filter(!dataset %in% c("brca","immuno","kirc","lihc","luad","prad","thca")) %>% 
  mutate(title = "Minimum absolute difference between groups for P <0.05",
         shp = dataset)


colnames(group_means.oth2)[1] <- "mean_value"
group_means.oth2$dataset <- case_when(group_means.oth2$dataset == "mts" ~ "Vaginal metatranscriptome",
                                     group_means.oth2$dataset == "sccyto" ~ "Single-cell transcriptome",
                                     group_means.oth2$dataset == "yeast" ~ "Yeast transcriptome")

group_means.yst <- group_means.oth2 %>% 
  filter(dataset == "Yeast transcriptome")

group_means.sc <- group_means.oth2 %>% 
  filter(dataset == "Single-cell transcriptome")

group_means.mts <- group_means.oth2 %>% 
  filter(dataset == "Vaginal metatranscriptome")

# plot stripchart: colours = tool
# png(paste0(repo,"figures/fig3_stripchartToolALDEx.png"),
#     units = "in", height = 4, width = 6, res = 600)

p2 <- ggplot(group_means.oth, aes(x = gamma.num, y = mean_value, colour=tool, shape = dataset))+
  geom_point(alpha = 1, size = 1)+
  geom_line(linewidth = 0.25)+
  scale_colour_manual(name = "Tool", values = c("dodgerblue2", "orangered2"))+
  stat_summary(data = group_means.hum, aes(group = tool), geom = "errorbar",
               fun.data = "mean_se", width = 0.0075, linewidth = 0.35, colour = "black")+
  stat_summary(data = group_means.hum, aes(group = tool, shape = shp),
               geom = "point", fun = "mean", size = 1.5)+
  stat_summary(data = group_means.hum, aes(group = tool),
               geom = "line", fun = "mean", linewidth = 0.25)+
  scale_x_continuous(limits = c(-0.005, 0.505), expand = c(0.0025,0.01))+
  scale_shape_manual(name = "Dataset",
                     values = c("Vaginal metatranscriptome" = 16, "Single-cell transcriptome" = 17,
                                "Yeast transcriptome" = 15, "TCGA transcriptomes" = 18))+
  labs(x = "Scale (\u03b3)", y = "Minimum log difference")+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.key.spacing.y = unit(0.1, "cm"), legend.key.spacing.x = unit(0, "cm"),
        legend.key.size = unit(0.1, "cm"), legend.margin = margin(0,0,0,0.0,"cm"),
        legend.title = element_text(face = "bold", size = 7), legend.text = element_text(size = 6),
        legend.box.spacing = unit(0.2, "cm"), strip.text = element_text(face = "bold"))

p2

# dev.off()


# plot stripchart: colours = dataset
# png(paste0(repo,"figures/fig3_stripchartDatasetALDEx.png"),
#     units = "in", height = 4, width = 6, res = 600)

p3 <- ggplot(group_means.oth2, aes(x = gamma.num, y = mean_value, colour=dataset, shape = tool))+
  labs(x = "Scale (\u03b3)", y = "Minimum log difference")+
  stat_summary(data = group_means.yst, aes(group = tool), geom = "errorbar",
               fun.data = "mean_se", width = 0.0055, linewidth = 0.35, colour = "black")+
  stat_summary(data = group_means.sc, aes(group = tool), geom = "errorbar",
               fun.data = "mean_se", width = 0.0055, linewidth = 0.35, colour = "black")+
  stat_summary(data = group_means.mts, aes(group = tool), geom = "errorbar",
               fun.data = "mean_se", width = 0.0055, linewidth = 0.35, colour = "black")+
  stat_summary(data = group_means.hum, aes(group = tool), geom = "errorbar",
               fun.data = "mean_se", width = 0.0055, linewidth = 0.35, colour = "black")+
  stat_summary(data = group_means.yst, aes(group = tool, colour = dataset),
               geom = "point", fun = "mean", size = 1)+
  stat_summary(data = group_means.sc, aes(group = tool, colour = dataset),
               geom = "point", fun = "mean", size = 1)+
  stat_summary(data = group_means.mts, aes(group = tool, colour = dataset),
               geom = "point", fun = "mean", size = 1)+
  stat_summary(data = group_means.hum, aes(group = tool,  colour = shp),
               geom = "point", fun = "mean", size = 1)+
  stat_summary(data = group_means.yst, aes(group = tool, colour = dataset),
               geom = "line", fun = "mean", linewidth = 0.25)+
  stat_summary(data = group_means.sc, aes(group = tool, colour = dataset),
               geom = "line", fun = "mean", linewidth = 0.25)+
  stat_summary(data = group_means.mts, aes(group = tool, colour = dataset),
               geom = "line", fun = "mean", linewidth = 0.25)+
  stat_summary(data = group_means.hum, aes(group = tool), colour = "purple3",
               geom = "line", fun = "mean", linewidth = 0.25)+
  scale_x_continuous(limits = c(-0.005, 0.505), expand = c(0.0025,0.01))+
  scale_colour_manual(name = "Dataset",
                      values = c("Vaginal metatranscriptome" = "dodgerblue", "Single-cell transcriptome" = "orangered2",
                                 "Yeast transcriptome" = "goldenrod2", "TCGA transcriptomes" = "purple3"))+
  labs(x = "Scale (\u03b3)", y = "Minimum log difference", shape = "Tool")+
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 9),
        legend.key.spacing.y = unit(0.1, "cm"), legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.size = unit(0.2, "cm"), legend.margin = margin(0,0,0,0.0,"cm"),
        legend.title = element_text(face = "bold", size = 8), legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.2, "cm"), strip.text = element_text(face = "bold"),
        legend.position = "top", legend.spacing.x = unit(0.5, "cm"))

p3

# dev.off()

# plot density plots and stripchart as a two-panel figure (strip = tool)
# png(paste0(repo,"figures/fig3_bothTool.png"),
#     units = "in", height = 4, width = 10, res = 600)

p1|p2

# dev.off()


# plot density plots and stripchart as a two-panel figure (strip = dataset)
# png(paste0(repo,"figures/fig3_bothDataset.png"),
#     units = "in", height = 4, width = 10, res = 600)

p1|p3

# dev.off()
