library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)


####new script 

anal.path <- "../ext_analysis/"
data <- "~/Documents/GitHub/ext_analysis/"
repo <- "~/Documents/GitHub/usri/"

source('code/summary.data.fun.R')

if(length(list.files(paste0(repo,"analysis/summarystats"), pattern = "ss")) != 88){
  
  # load in analysis results (takes a few minutes!)
  for(i in list.files(data, pattern = "immuno.*aldex")){
    load(paste0(data, i))
  }
  
  # build summary object from datasets (will take several minutes)
  i <- ls(pattern = "immuno")[1]
    
  # extract dataset and tool names from input file (sensitive to aldex gamma)
  dataset <- str_split(i, "\\." , 3)[[1]][1]
  gamma <- "data_"
  tool <- "aldex"
    
  # run 'get_confusion.R' on analysis objects to automatically build the
  # confusion matrices
  assign(x = paste0("ss.", dataset),
         sum.fun(aldex2_0 = get(paste0(dataset, ".", gamma, 0, ".", tool, 2)),
                 aldex2_1 = get(paste0(dataset, ".", gamma, 1, ".", tool, 2)),
                 aldex2_2 = get(paste0(dataset, ".", gamma, 2, ".", tool, 2)),
                 aldex2_5 = get(paste0(dataset, ".", gamma, 5, ".", tool, 2)),
                 aldex3_0 = get(paste0(dataset, ".", gamma, 0, ".", tool, 3)),
                 aldex3_1 = get(paste0(dataset, ".", gamma, 1, ".", tool, 3)),
                 aldex3_2 = get(paste0(dataset, ".", gamma, 2, ".", tool, 3)),
                 aldex3_5 = get(paste0(dataset, ".", gamma, 5, ".", tool, 3))))
  
  
  # save confusion matrix list objects as .Rda
  
    save(list = paste0("ss.", dataset),
         file = paste0(repo,"analysis/summarystats/ss.", dataset, ".Rda"))
  
  # extract TPFPR matrices from each confusion matrix and delete 'cm.' objects
  for(i in ls(pattern = "ss")){
    assign(x = paste0("minmax", gsub("cm\\.", "", i)),
           value = as.data.frame(get(i)$TPFPR))
    rm(list = i)
  }
  
  # finally, remove analysis objects from environment
  for(i in ls(pattern = "^immuno")){
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

# # filter dataframe for TPR & FDR values (using the different definitions of
# # true negatives and with/without a threshold of 0.5)
# plot.cFDR0 <- plot.df %>% 
#   filter(metric == "cFDR0")
# 
# plot.cTPR0 <- plot.df %>% 
#   filter(metric == "cTPR0")
# 
# plot.cFDR5 <- plot.df %>% 
#   filter(metric == "cFDR5")
# 
# plot.cTPR5 <- plot.df %>% 
#   filter(metric == "cTPR5")
# 
# plot.zFDR0 <- plot.df %>% 
#   filter(metric == "zFDR0")
# 
# plot.zTPR0 <- plot.df %>% 
#   filter(metric == "zTPR0")
# 
# plot.zFDR5 <- plot.df %>% 
#   filter(metric == "zFDR5")
# 
# plot.zTPR5 <- plot.df %>% 
#   filter(metric == "zTPR5")

############################## plotting function ##############################

# define function for plotting data: all tools on same graph and faceted by 
# dataset. User defines which metric is to be used for plotting (i.e. plot TPR
# or FDR with which definition of true negatives (coefficient vs. zero) and
# with/without an arbitrary coefficient threshold of 0.5 for positives)

function(df = NULL, plot.type = NULL, tn.def = NULL, coeff.thresh = NULL){
  
  # throw errors if inputs are not correct
  if(is.null(df)) stop("Please specify dataframe containing plotting data")
  if(is.null(plot.type)) stop("Please specify metric to plot: TPR or FDR")
  if(is.null(tn.def)) stop("Please specify which definition of true negatives to use: 'coeff' or 'zero'")
  if(is.null(coeff.thresh)) stop("Please specify if a coefficient threshold of 0.5 should be used for defining positives")
  
  if(!is.data.frame(df)) stop("Argument 'df' must be a dataframe")
  
  if(!plot.type %in% c("FDR", "TPR")) stop("Argument 'plot.type' must be one of: FDR, TPR")
  if(length(plot.type) != 1) stop("Argument 'plot.type' must be a vector of length = 1")
  if(!is.vector(plot.type) != 1) stop("Argument 'plot.type' must be a vector of length = 1")
  
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
  
  # plot data
  to.plot %>% 
    ggplot(aes(x = coeff, y = value))
}


  
  
  
  
  
  
  
  
  
  
load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

imm.sum <- sum.fun(aldex2_0 = immuno.data_0.aldex2, aldex2_1 = immuno.data_1.aldex2, aldex2_2 = immuno.data_2.aldex2, aldex2_5 = immuno.data_5.aldex2, aldex3_0 = immuno.data_0.aldex3, aldex3_1 = immuno.data_1.aldex3, aldex3_2 = immuno.data_2.aldex3, aldex3_5 = immuno.data_5.aldex3)

##########
#Immuno Plot -----------------------------------------
df1 <- as.data.frame(imm.sum$minmax2_0)
df2 <- as.data.frame(imm.sum$minmax2_1)
df3 <- as.data.frame(imm.sum$minmax2_2)
df4 <- as.data.frame(imm.sum$minmax2_5)
df5 <- as.data.frame(imm.sum$minmax3_0)
df6 <- as.data.frame(imm.sum$minmax3_1)
df7 <- as.data.frame(imm.sum$minmax3_2)
df8 <- as.data.frame(imm.sum$minmax3_5)
df1$scale<-'γ = 0'
df2$scale<-'γ = 0.1'
df3$scale<-'γ = 0.2'
df4$scale<-'γ = 0.5'
df5$scale<-'γ = 0'
df6$scale<-'γ = 0.1'
df7$scale<-'γ = 0.2'
df8$scale<-'γ = 0.5'
df1$tool<-'ALDEx2'
df2$tool<-'ALDEx2'
df3$tool<-'ALDEx2'
df4$tool<-'ALDEx2'
df5$tool<-'ALDEx3'
df6$tool<-'ALDEx3'
df7$tool<-'ALDEx3'
df8$tool<-'ALDEx3'
df_merged = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8)

group_means <- df_merged %>%
  group_by(tool, scale) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

df_merged %>% 
  ggplot(aes(x = abs, fill=scale)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=scale), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')

######################

load(paste(anal.path,"prad.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_5.Rda", sep=""))

prad.sum <- sum.fun(aldex2_0 = prad.data_0.aldex2, aldex2_1 = prad.data_1.aldex2, aldex2_2 = prad.data_2.aldex2, aldex2_5 = prad.data_5.aldex2, aldex3_0 = prad.data_0.aldex3, aldex3_1 = prad.data_1.aldex3, aldex3_2 = prad.data_2.aldex3, aldex3_5 = prad.data_5.aldex3)

#PRAD plot ------------------------------


df1 <- as.data.frame(prad.sum$minmax2_0)
df2 <- as.data.frame(prad.sum$minmax2_1)
df3 <- as.data.frame(prad.sum$minmax2_2)
df4 <- as.data.frame(prad.sum$minmax2_5)
df5 <- as.data.frame(prad.sum$minmax3_0)
df6 <- as.data.frame(prad.sum$minmax3_1)
df7 <- as.data.frame(prad.sum$minmax3_2)
df8 <- as.data.frame(prad.sum$minmax3_5)
df1$scale<-'γ = 0'
df2$scale<-'γ = 0.1'
df3$scale<-'γ = 0.2'
df4$scale<-'γ = 0.5'
df5$scale<-'γ = 0'
df6$scale<-'γ = 0.1'
df7$scale<-'γ = 0.2'
df8$scale<-'γ = 0.5'
df1$tool<-'ALDEx2'
df2$tool<-'ALDEx2'
df3$tool<-'ALDEx2'
df4$tool<-'ALDEx2'
df5$tool<-'ALDEx3'
df6$tool<-'ALDEx3'
df7$tool<-'ALDEx3'
df8$tool<-'ALDEx3'
df_merged = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8)

group_means <- df_merged %>%
  group_by(tool, scale) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

df_merged %>% 
  ggplot(aes(x = abs, fill=scale)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=scale), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')

######################



######################

load(paste(anal.path,"luad.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_5.Rda", sep=""))

luad.sum <- sum.fun(aldex2_0 = luad.data_0.aldex2, aldex2_1 = luad.data_1.aldex2, aldex2_2 = luad.data_2.aldex2, aldex2_5 = luad.data_5.aldex2, aldex3_0 = luad.data_0.aldex3, aldex3_1 = luad.data_1.aldex3, aldex3_2 = luad.data_2.aldex3, aldex3_5 = luad.data_5.aldex3)

#PRAD plot ------------------------------
df1 <- as.data.frame(luad.sum$minmax2_0)
df2 <- as.data.frame(luad.sum$minmax2_1)
df3 <- as.data.frame(luad.sum$minmax2_2)
df4 <- as.data.frame(luad.sum$minmax2_5)
df5 <- as.data.frame(luad.sum$minmax3_0)
df6 <- as.data.frame(luad.sum$minmax3_1)
df7 <- as.data.frame(luad.sum$minmax3_2)
df8 <- as.data.frame(luad.sum$minmax3_5)
df1$scale<-'γ = 0'
df2$scale<-'γ = 0.1'
df3$scale<-'γ = 0.2'
df4$scale<-'γ = 0.5'
df5$scale<-'γ = 0'
df6$scale<-'γ = 0.1'
df7$scale<-'γ = 0.2'
df8$scale<-'γ = 0.5'
df1$tool<-'ALDEx2'
df2$tool<-'ALDEx2'
df3$tool<-'ALDEx2'
df4$tool<-'ALDEx2'
df5$tool<-'ALDEx3'
df6$tool<-'ALDEx3'
df7$tool<-'ALDEx3'
df8$tool<-'ALDEx3'
df_merged = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8)

group_means <- df_merged %>%
  group_by(tool, scale) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

df_merged %>% 
  ggplot(aes(x = abs, fill=scale)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=scale), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')

######################


######################

load(paste(anal.path,"thca.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_5.Rda", sep=""))

thca.sum <- sum.fun(aldex2_0 = thca.data_0.aldex2, aldex2_1 = thca.data_1.aldex2, aldex2_2 = thca.data_2.aldex2, aldex2_5 = thca.data_5.aldex2, aldex3_0 = thca.data_0.aldex3, aldex3_1 = thca.data_1.aldex3, aldex3_2 = thca.data_2.aldex3, aldex3_5 = thca.data_5.aldex3)

#PRAD plot ------------------------------
df1 <- as.data.frame(thca.sum$minmax2_0)
df2 <- as.data.frame(thca.sum$minmax2_1)
df3 <- as.data.frame(thca.sum$minmax2_2)
df4 <- as.data.frame(thca.sum$minmax2_5)
df5 <- as.data.frame(thca.sum$minmax3_0)
df6 <- as.data.frame(thca.sum$minmax3_1)
df7 <- as.data.frame(thca.sum$minmax3_2)
df8 <- as.data.frame(thca.sum$minmax3_5)
df1$scale<-'γ = 0'
df2$scale<-'γ = 0.1'
df3$scale<-'γ = 0.2'
df4$scale<-'γ = 0.5'
df5$scale<-'γ = 0'
df6$scale<-'γ = 0.1'
df7$scale<-'γ = 0.2'
df8$scale<-'γ = 0.5'
df1$tool<-'ALDEx2'
df2$tool<-'ALDEx2'
df3$tool<-'ALDEx2'
df4$tool<-'ALDEx2'
df5$tool<-'ALDEx3'
df6$tool<-'ALDEx3'
df7$tool<-'ALDEx3'
df8$tool<-'ALDEx3'
df_merged = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8)

group_means <- df_merged %>%
  group_by(tool, scale) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

df_merged %>% 
  ggplot(aes(x = abs, fill=scale)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=scale), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')

######################


######################

load(paste(anal.path,"lihc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_5.Rda", sep=""))

lihc.sum <- sum.fun(aldex2_0 = lihc.data_0.aldex2, aldex2_1 = lihc.data_1.aldex2, aldex2_2 = lihc.data_2.aldex2, aldex2_5 = lihc.data_5.aldex2, aldex3_0 = lihc.data_0.aldex3, aldex3_1 = lihc.data_1.aldex3, aldex3_2 = lihc.data_2.aldex3, aldex3_5 = lihc.data_5.aldex3)

#PRAD plot ------------------------------
df1 <- as.data.frame(lihc.sum$minmax2_0)
df2 <- as.data.frame(lihc.sum$minmax2_1)
df3 <- as.data.frame(lihc.sum$minmax2_2)
df4 <- as.data.frame(lihc.sum$minmax2_5)
df5 <- as.data.frame(lihc.sum$minmax3_0)
df6 <- as.data.frame(lihc.sum$minmax3_1)
df7 <- as.data.frame(lihc.sum$minmax3_2)
df8 <- as.data.frame(lihc.sum$minmax3_5)
df1$scale<-'γ = 0'
df2$scale<-'γ = 0.1'
df3$scale<-'γ = 0.2'
df4$scale<-'γ = 0.5'
df5$scale<-'γ = 0'
df6$scale<-'γ = 0.1'
df7$scale<-'γ = 0.2'
df8$scale<-'γ = 0.5'
df1$tool<-'ALDEx2'
df2$tool<-'ALDEx2'
df3$tool<-'ALDEx2'
df4$tool<-'ALDEx2'
df5$tool<-'ALDEx3'
df6$tool<-'ALDEx3'
df7$tool<-'ALDEx3'
df8$tool<-'ALDEx3'
df_merged = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8)

group_means <- df_merged %>%
  group_by(tool, scale) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

df_merged %>% 
  ggplot(aes(x = abs, fill=scale)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=scale), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')

######################
#Plot with stripchart

##########
#Immuno Plot -----------------------------------------
df1 <- as.data.frame(imm.sum$minmax2_0)
df2 <- as.data.frame(imm.sum$minmax2_1)
df3 <- as.data.frame(imm.sum$minmax2_2)
df4 <- as.data.frame(imm.sum$minmax2_5)
df5 <- as.data.frame(imm.sum$minmax3_0)
df6 <- as.data.frame(imm.sum$minmax3_1)
df7 <- as.data.frame(imm.sum$minmax3_2)
df8 <- as.data.frame(imm.sum$minmax3_5)
df1$scale<-c(0)
df2$scale<-c(0.1)
df3$scale<-c(0.2)
df4$scale<-c(0.5)
df5$scale<-c(0)
df6$scale<-c(0.1)
df7$scale<-c(0.2)
df8$scale<-c(0.5)
df1$scale1<-'0'
df2$scale1<-'0.1'
df3$scale1<-'0.2'
df4$scale1<-'0.5'
df5$scale1<-'0'
df6$scale1<-'0.1'
df7$scale1<-'0.2'
df8$scale1<-'0.5'
df1$tool<-'ALDEx2'
df2$tool<-'ALDEx2'
df3$tool<-'ALDEx2'
df4$tool<-'ALDEx2'
df5$tool<-'ALDEx3'
df6$tool<-'ALDEx3'
df7$tool<-'ALDEx3'
df8$tool<-'ALDEx3'
df1$Dataset<-'Immuno'
df2$Dataset<-'Immuno'
df3$Dataset<-'Immuno'
df4$Dataset<-'Immuno'
df5$Dataset<-'Immuno'
df6$Dataset<-'Immuno'
df7$Dataset<-'Immuno'
df8$Dataset<-'Immuno'
df1l <- as.data.frame(lihc.sum$minmax2_0)
df2l <- as.data.frame(lihc.sum$minmax2_1)
df3l <- as.data.frame(lihc.sum$minmax2_2)
df4l <- as.data.frame(lihc.sum$minmax2_5)
df5l <- as.data.frame(lihc.sum$minmax3_0)
df6l <- as.data.frame(lihc.sum$minmax3_1)
df7l <- as.data.frame(lihc.sum$minmax3_2)
df8l <- as.data.frame(lihc.sum$minmax3_5)
df1l$scale<-c(0)
df2l$scale<-c(0.1)
df3l$scale<-c(0.2)
df4l$scale<-c(0.5)
df5l$scale<-c(0)
df6l$scale<-c(0.1)
df7l$scale<-c(0.2)
df8l$scale<-c(0.5)
df1l$scale1<-'0'
df2l$scale1<-'0.1'
df3l$scale1<-'0.2'
df4l$scale1<-'0.5'
df5l$scale1<-'0'
df6l$scale1<-'0.1'
df7l$scale1<-'0.2'
df8l$scale1<-'0.5'
df1l$tool<-'ALDEx2'
df2l$tool<-'ALDEx2'
df3l$tool<-'ALDEx2'
df4l$tool<-'ALDEx2'
df5l$tool<-'ALDEx3'
df6l$tool<-'ALDEx3'
df7l$tool<-'ALDEx3'
df8l$tool<-'ALDEx3'
df1l$Dataset<-'LIHC'
df2l$Dataset<-'LIHC'
df3l$Dataset<-'LIHC'
df4l$Dataset<-'LIHC'
df5l$Dataset<-'LIHC'
df6l$Dataset<-'LIHC'
df7l$Dataset<-'LIHC'
df8l$Dataset<-'LIHC'
df1p <- as.data.frame(prad.sum$minmax2_0)
df2p <- as.data.frame(prad.sum$minmax2_1)
df3p <- as.data.frame(prad.sum$minmax2_2)
df4p <- as.data.frame(prad.sum$minmax2_5)
df5p <- as.data.frame(prad.sum$minmax3_0)
df6p <- as.data.frame(prad.sum$minmax3_1)
df7p <- as.data.frame(prad.sum$minmax3_2)
df8p <- as.data.frame(prad.sum$minmax3_5)
df1p$scale<-c(0)
df2p$scale<-c(0.1)
df3p$scale<-c(0.2)
df4p$scale<-c(0.5)
df5p$scale<-c(0)
df6p$scale<-c(0.1)
df7p$scale<-c(0.2)
df8p$scale<-c(0.5)
df1p$scale1<-'0'
df2p$scale1<-'0.1'
df3p$scale1<-'0.2'
df4p$scale1<-'0.5'
df5p$scale1<-'0'
df6p$scale1<-'0.1'
df7p$scale1<-'0.2'
df8p$scale1<-'0.5'
df1p$tool<-'ALDEx2'
df2p$tool<-'ALDEx2'
df3p$tool<-'ALDEx2'
df4p$tool<-'ALDEx2'
df5p$tool<-'ALDEx3'
df6p$tool<-'ALDEx3'
df7p$tool<-'ALDEx3'
df8p$tool<-'ALDEx3'
df1p$Dataset<-'PRAD'
df2p$Dataset<-'PRAD'
df3p$Dataset<-'PRAD'
df4p$Dataset<-'PRAD'
df5p$Dataset<-'PRAD'
df6p$Dataset<-'PRAD'
df7p$Dataset<-'PRAD'
df8p$Dataset<-'PRAD'

df1t <- as.data.frame(thca.sum$minmax2_0)
df2t <- as.data.frame(thca.sum$minmax2_1)
df3t <- as.data.frame(thca.sum$minmax2_2)
df4t <- as.data.frame(thca.sum$minmax2_5)
df5t <- as.data.frame(thca.sum$minmax3_0)
df6t <- as.data.frame(thca.sum$minmax3_1)
df7t <- as.data.frame(thca.sum$minmax3_2)
df8t <- as.data.frame(thca.sum$minmax3_5)
df1t$scale<-c(0)
df2t$scale<-c(0.1)
df3t$scale<-c(0.2)
df4t$scale<-c(0.5)
df5t$scale<-c(0)
df6t$scale<-c(0.1)
df7t$scale<-c(0.2)
df8t$scale<-c(0.5)
df1t$scale1<-'0'
df2t$scale1<-'0.1'
df3t$scale1<-'0.2'
df4t$scale1<-'0.5'
df5t$scale1<-'0'
df6t$scale1<-'0.1'
df7t$scale1<-'0.2'
df8t$scale1<-'0.5'
df1t$tool<-'ALDEx2'
df2t$tool<-'ALDEx2'
df3t$tool<-'ALDEx2'
df4t$tool<-'ALDEx2'
df5t$tool<-'ALDEx3'
df6t$tool<-'ALDEx3'
df7t$tool<-'ALDEx3'
df8t$tool<-'ALDEx3'
df1t$Dataset<-'THCA'
df2t$Dataset<-'THCA'
df3t$Dataset<-'THCA'
df4t$Dataset<-'THCA'
df5t$Dataset<-'THCA'
df6t$Dataset<-'THCA'
df7t$Dataset<-'THCA'
df8t$Dataset<-'THCA'

df_merged = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8)
df_merged0 = bind_rows(df1,df2,df3,df4, df5, df6, df7, df8, df1l, df2l, df3l, df4l, df5l, df6l, df7l, df8l, df1p,df2p,df3p,df4p, df5p, df6p, df7p, df8p, df1t, df2t, df3t, df4t, df5t, df6t, df7t, df8t)

group_means <- df_merged %>%
  group_by(tool, scale1) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

group_means0 <- df_merged0 %>%
  group_by(tool, scale, Dataset) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

p1<-df_merged %>% 
  ggplot(aes(x = abs, fill=scale1)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=scale1), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  labs(title = "Immunotherapy Minimum Difference Density Distribution")+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')


p2 <- ggplot(group_means0, aes(x = scale, y = mean_value, colour=Dataset)) +
  geom_jitter(width = 0.01, alpha = 0.7) + 
  labs(title = "Average Minimum Difference Across Datasets",
       x = "Scale",
       y = "Mean Absolute Minimum Difference for Significance") +
  theme_bw()

p1|p2
