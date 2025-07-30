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

# NOTE: as long as all 132 confusion matrix objects are present within the 
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

if(length(list.files(paste0(repo,"analysis/confusionMats"), pattern = "cm")) != 77){
  
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
                         FDR = 0.1)) 
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

# collapse list to df and delete list
plot.df <- do.call(rbind, tmplist)
rm(tmplist)








