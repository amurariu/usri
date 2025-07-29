# summarising FDRs and TPRs across all tools & datasets

# Scott Dos Santos
# Last edited: 2025-07-29

# script to build summary dataframes and graphs for showing average TPR and FDR
# of all tools at varying prescribed coefficient thresholds, per dataset. Uses
# pre-computed TPR and FDR calculated from the confusion matrix script in this
# repo ('get_confusion.R').

#################################### setup ####################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# local directory where analysis results live (all too large for GH)
data <- "~/Documents/GitHub/ext_analysis/"

# path to repo directory
repo <- "~/Documents/GitHub/usri/"

# source confusion matrix code
source(paste0(repo, "code/get_confusion.R"))

############################## confusion matrices ##############################

# load in analysis data, get confusion matrices for all objects and extract the 
# TPFPR matrix summary (or load them if they already exist as .Rda files) 

if(length(list.files(paste0(repo,"analysis/confusionMats"), pattern = "cm")) != 11){
  
  # load in analysis results (takes a minute or so!)
  for(i in list.files(data, pattern = "immuno")){
    load(paste0(data, i))
  }
  
  # build confusion matrix from datasets
  cm.immuno.aldex2.0 <- get_confusion(input = immuno.data_0.aldex2, prog = "ALDEx2", FDR = 0.1)
  cm.immuno.aldex2.1 <- get_confusion(input = immuno.data_1.aldex2, prog = "ALDEx2", FDR = 0.1)
  cm.immuno.aldex2.2 <- get_confusion(input = immuno.data_2.aldex2, prog = "ALDEx2", FDR = 0.1)
  cm.immuno.aldex2.5 <- get_confusion(input = immuno.data_5.aldex2, prog = "ALDEx2", FDR = 0.1)
  cm.immuno.aldex3.0 <- get_confusion(input = immuno.data_0.aldex3, prog = "ALDEx3", FDR = 0.1)
  cm.immuno.aldex3.1 <- get_confusion(input = immuno.data_1.aldex3, prog = "ALDEx3", FDR = 0.1)
  cm.immuno.aldex3.2 <- get_confusion(input = immuno.data_2.aldex3, prog = "ALDEx3", FDR = 0.1)
  cm.immuno.aldex3.5 <- get_confusion(input = immuno.data_5.aldex3, prog = "ALDEx3", FDR = 0.1)
  cm.immuno.deseq <- get_confusion(input = immuno.data.DESeq, prog = "DESeq", FDR = 0.1)
  cm.immuno.edger <- get_confusion(input = immuno.data.edgeR, prog = "edgeR", FDR = 0.1)
  cm.immuno.limma <- get_confusion(input = immuno.data.limma, prog = "limma", FDR = 0.1)
  
  # save confusion matrix list objects
  for(i in ls(pattern = "cm")){
    save(list = i,
         file = paste0(repo,"analysis/confusionMats/", i, ".Rda"))
  }
  
  # extract TPFPR matrices from each confusion matrix
  for(i in ls(pattern = "cm")){
    assign(x = paste0("tpfpr.", gsub("cm\\.", "", i)),
           value = as.data.frame(get(i)$TPFPR))
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

# collapse list to df and delete list
plot.df <- do.call(rbind, tmplist)
rm(tmplist)








