# summarising FDRs and TPRs across all tools & datasets

# Scott Dos Santos
# Last edited: 2025-08-08

# script to build summary dataframes and graphs for showing average TPR and FDR
# of all tools at varying prescribed coefficient thresholds, per dataset. Uses
# pre-computed TPR and FDR calculated from the confusion matrix script in this
# repo ('get_confusion.R').

#################################### setup ####################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

# local directory where analysis results live (all too large for GH)
data <- "~/Documents/GitHub/ext_analysis/"

# path to repo directory
repo <- "~/Documents/GitHub/usri/"

# source confusion matrix code
source(paste0(repo, "code/get_confusion.R"))

############################## confusion matrices ##############################

# load in analysis data, get confusion matrices for all objects and extract the 
# TPFPR matrix summary (or load them if they already exist as .Rda files) 

# NOTE: as long as all confusion matrix objects are present within the 
#       'analysis/confusionMats/' directory, this loop should just load all of
#       the 'TPFPR' matrices for each dataset and tool. To get around memory
#       limits on my local machine, I had to run the code in the 'if' loop below
#       manually for each dataset (changing the three 'pattern = ' arguments)
#       to generate the confusion matrix objects. This shouldn't be a problem
#       for anyone looking at this file after me who wants to replicate the
#       analysis we did (as these files are already generated). If for whatever
#       reason you want to re-make the confusion matrix objects, you'll have to
#       re-run the if statement code manually. Sorry (not sorry) for being lazy.
#           - Scott Dos Santos, 30th July 2025.

if(length(list.files(paste0(repo,"analysis/confusionMats"), pattern = "cm")) != 150){
  
  # load in analysis results (takes a few minutes!)
  for(i in list.files(data, pattern = "thca")){
    load(paste0(data, i))
  }
  
  # intialise vectors for storing dataset and tool strings
  dataset <- vector()
  tool <- vector()
  
  # build confusion matrix from datasets (will take several minutes)
  for(i in ls(pattern = "thca")){
    
    # extract dataset and tool names from input file (sensitive to aldex gamma)
    dataset <- str_split(i, "\\." , 3)[[1]][1]
    tool <- tolower(str_split(i, "\\." , 3)[[1]][3])
    if(tool == "aldex2"| tool == "aldex3"){
      tool <- paste0(tool, ".", gsub("data_", "", str_split(i, "\\." , 3)[[1]][2]))
    }
    
    # run 'get_confusion.R' on analysis objects to automatically build the
    # confusion matrices
    assign(x = paste0("cm.", dataset, ".", tool),
           get_confusion(input = get(i),
                         prog = gsub("\\..", "", tool),
                         FDR = 0.05)) 
  }
  
  # save confusion matrix list objects as .Rda
  for(i in ls(pattern = "cm")){
    save(list = i,
         file = paste0(repo,"analysis/confusionMats/", i, ".Rda"))
  }
  
  # extract TPFPR matrices from each confusion matrix and delete 'cm.' objects
  for(i in ls(pattern = "cm")){
    assign(x = paste0("tpfpr.", gsub("cm\\.", "", i)),
           value = as.data.frame(get(i)$TPFPR))
    rm(list = i)
  }
  
  # finally, remove analysis objects from environment
  for(i in ls(pattern = "^thca")){
    rm(list = i)
  }
} else{
  
  # load in confusion matrix objects from .Rda
  for(i in list.files(paste0(repo,"analysis/confusionMats"), pattern = "cm")){
    load(paste0(repo, "analysis/confusionMats/", i))
  }
  
  # extract TPFPR matrices from each confusion matrix
  for(i in ls(pattern = "cm")){
    assign(x = paste0("tpfpr.", gsub("cm\\.", "", i)),
           value = as.data.frame(get(i)$TPFPR))
    rm(list = i)
  }
}

############################# data transformation #############################

# transform all data in tpfpr to tidy format for ggplot
tmplist <- list()
for(i in ls(pattern = "tpfpr")){
  # extract rownames and set to 'coeff'
  df <- get(i)
  df$coeff <- rownames(df)
  rownames(df) <- NULL
  
  # pivot around coeff column and add columns indicating dataset and tool
  df <- df %>% 
    pivot_longer(-coeff, names_to = "metric", values_to = "value") %>% 
    arrange(metric, coeff)
  
  df$dataset <- str_split(i, "\\.", 3)[[1]][2]
  df$tool <- str_split(i, "\\.", 3)[[1]][3]
  
  # add to list
  tmplist[[i]] <-df
  rm(df)
}

# collapse list to df, then remove rownames, convert 'coeff' to numeric and
# change order of columns
plot.df <- do.call(rbind, tmplist)
plot.df$coeff <- as.numeric(plot.df$coeff)
plot.df <- plot.df %>% 
  relocate(c(dataset, tool), .before = coeff) %>% 
  relocate(metric, .after = tool)

# delete tpfpr objects and temporary objects
rm(list = c("tmplist", "i", ls(pattern = "tpfpr")))

# edit tool column to have sentence case and show gamma symbols
plot.df$tool <- case_when(plot.df$tool == "deseq" ~ "DESeq2",
                          plot.df$tool == "edger" ~ "EdgeR",
                          plot.df$tool == "limma" ~ "Limma",
                          .default = plot.df$tool)

plot.df$tool <- gsub("aldex", "ALDEx", plot.df$tool)
plot.df$tool <- gsub("\\.0", " (\u03b3 = 0)", plot.df$tool)
plot.df$tool <- gsub("\\.1", " (\u03b3 = 0.1)", plot.df$tool)
plot.df$tool <- gsub("\\.2", " (\u03b3 = 0.2)", plot.df$tool)
plot.df$tool <- gsub("\\.3", " (\u03b3 = 0.3)", plot.df$tool)
plot.df$tool <- gsub("\\.4", " (\u03b3 = 0.4)", plot.df$tool)
plot.df$tool <- gsub("\\.5", " (\u03b3 = 0.5)", plot.df$tool)

# edit dataset column for sentence case
plot.df$dataset <- case_when(plot.df$dataset == "brca" ~ "Cancer genome atlas: BRCA",
                             plot.df$dataset == "immuno" ~ "Immunotherapy transcriptome",
                             plot.df$dataset == "kirc" ~ "Cancer genome atlas: KIRC",
                             plot.df$dataset == "lihc" ~ "Cancer genome atlas: LIHC",
                             plot.df$dataset == "luad" ~ "Cancer genome atlas: LUAD",
                             plot.df$dataset == "mts" ~ "Vaginal metatranscriptome",
                             plot.df$dataset == "prad" ~ "Cancer genome atlas: PRAD",
                             plot.df$dataset == "sccyto" ~ "Single-cell transcriptome",
                             plot.df$dataset == "thca" ~ "Cancer genome atlas: THCA",
                             plot.df$dataset == "yeast" ~ "Yeast transcriptome",
                             .default = plot.df$dataset)

# convert tool column to factor and change levels so that aldex2/3 gamma values
# are in order from 0 to 0.5
lvl <- levels(factor(plot.df$tool))

plot.df$tool <- factor(plot.df$tool, 
                       levels = c(lvl[6],lvl[1:5], lvl[12], lvl[7:11], lvl[13:15]))

############################## plotting function ##############################

# define function for plotting data: all tools on same graph and faceted by 
# dataset. User defines which metric is to be used for plotting (i.e. plot TPR
# or FDR with which definition of true negatives (coefficient vs. zero) and
# with/without an arbitrary coefficient threshold of 0.5 for positives)

nicePlots <- function(df = NULL, plot.type = NULL, tn.def = NULL, coeff.thresh = NULL){
  
  # throw errors if inputs are not correct
  if(is.null(df)) stop("Please specify dataframe containing plotting data")
  if(is.null(plot.type)) stop("Please specify metric to plot: TPR or FDR")
  if(is.null(tn.def)) stop("Please specify which definition of true negatives to use: 'coeff' or 'zero'")
  if(is.null(coeff.thresh)) stop("Please specify if a coefficient threshold of 0.5 should be used for defining positives")
  
  if(!is.data.frame(df)) stop("Argument 'df' must be a dataframe")
  
  if(!plot.type %in% c("FDR", "TPR")) stop("Argument 'plot.type' must be one of: FDR, TPR")
  if(length(plot.type) != 1) stop("Argument 'plot.type' must be a vector of length = 1")
  if(!is.vector(plot.type)) stop("Argument 'plot.type' must be a vector of length = 1")
  
  if(!tn.def %in% c("coeff", "zero")) stop("Argument 'tn.def' must be one of: coeff, zero")
  if(length(tn.def) != 1) stop("Argument 'tn.def' must be a vector of length = 1")
  if(!is.vector(tn.def)) stop("Argument 'tn.def' must be a vector of length = 1")
  
  if(!coeff.thresh %in% c(TRUE, FALSE)) stop("Argument 'coeff.thresh' must be TRUE or FALSE")
  if(length(coeff.thresh) != 1) stop("Argument 'coeff.thresh' must be logical and of length = 1")
  if(!is.logical(coeff.thresh)) stop("Argument 'coeff.thresh' must be logical")
  
  # filter plotting data based on inputs to function
  plot.a <- plot.type
  plot.b <- switch(tn.def, "coeff" = "c", "zero" = "z")
  plot.c <- switch(as.character(coeff.thresh), "TRUE" = "5", "FALSE" = "0")
  to.filter <- paste0(plot.b, plot.a, plot.c)
  
  to.plot <- df %>% 
    filter(metric == to.filter)
  
  # get y axis label and title based on plot type input
  plot.ylab <- switch(plot.type,
                      "FDR" = "FDR",
                      "TPR" = "TPR")
  
  # get y axis label and title based on plot type input
  plot.title <- switch(plot.type,
                       "FDR" = "False discovery rate (FDR) by dataset and tool",
                       "TPR" = "Sensitivity / true positive rate (TPR) by dataset and tool")
  
  # define colour scheme for datasets
  cols.ald2 <- c("#E6E6FA","#CBE3FC","#B0E2FF","#66B3F6","#2171B5","#08306B")
  cols.ald3 <- brewer.pal(n = 9, name = "YlOrRd")[c(1,3,5,6,8,9)]
  
  cols.dataset <- c(cols.ald2,    # ALDEx2, gamma = 0 - 0.5
                    cols.ald3,    # ALDEx3, gamma = 0 - 0.5
                    "black",      # DESeq2
                    "grey50",     # EdgeR
                    "grey80"      # limma
  )
  
  # plot data
  nice.plot <- to.plot %>% 
    ggplot(aes(x = coeff, y = value, colour = tool))+
    geom_point(size = 1)+
    geom_line()+
    xlab("Threshold: model difference between groups")+
    ylab(plot.ylab)+
    # ggtitle(plot.title)+
    scale_colour_manual(name = "Tool", values = cols.dataset)+
    theme_bw()+
    guides(colour = guide_legend(nrow = 2))+
    theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10),
          legend.key.spacing.y = unit(0.1, "cm"), legend.key.spacing.x = unit(0, "cm"),
          legend.key.size = unit(0.3, "cm"), legend.margin = margin(0,0,0,0.0,"cm"),
          legend.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 9),
          legend.box.spacing = unit(0.2, "cm"), strip.text = element_text(face = "bold"),
          legend.position = "top")+
    facet_wrap(~dataset, ncol = 5)
  
  return(nice.plot)
}

############################## plotting TPR & FDR ##############################

# NOTE: these figures need to be edited in Inkscape to change the x axis 0 value
#       from '0.00' to '0', to avoid spacing issues

# false discovery rate: negatives >coeff, no difference threshold (cFDR0)
# png(paste0(repo, "figures/fig1_tprfdr_falseDiscoveryRate.png"),
#     units = "in", height = 6, width = 12, res = 600)

nicePlots(df = plot.df, plot.type = "FDR", tn.def = "coeff", coeff.thresh = FALSE)

# dev.off()

# sensitivity: negatives >coeff, no difference threshold (cTPR0)
# png(paste0(repo, "figures/fig1_tprfdr_sensitivity.png"),
#     units = "in", height = 6, width = 12, res = 600)

nicePlots(df = plot.df, plot.type = "TPR", tn.def = "coeff", coeff.thresh = FALSE)

# dev.off()
