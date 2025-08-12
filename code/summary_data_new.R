library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(ALDEx2,warn.conflicts = F) 


anal.path <- "../ext_analysis/"
data <- "~/Documents/GitHub/ext_analysis/"
repo <- "~/Documents/GitHub/usri/"

source('code/summary.data.fun.R')

if(length(list.files(paste0(repo,"analysis/summarystats"), pattern = "ss")) != 7){
  
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
  
  # finally, remove analysis objects from environment
  for(i in ls(pattern = "^immuno")){
    rm(list = i)
  }
} else{
  
  # load in confusion matrix objects from .Rda
  for(i in list.files(paste0(repo,"analysis/summarystats"), pattern = "ss")){
    load(paste0(repo, "analysis/summarystats/", i))
  }
}

############################# data transformation #############################

tmplist <- list()

for(i in ls(pattern = "^ss.")){
  # pull list
  df <- get(i)
  
  # make a temp dataframe with 1600 rows (200 x 8 tool/gamma combos) and 3 cols:
  # abs values, dataset and tool
  tempdf <- data.frame(matrix(data = NA, nrow = 1600, ncol = 3))
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
}

# collapse list to df, then remove rownames, convert 'coeff' to numeric and
# change order of columns
plot.df <- do.call(rbind, tmplist)

# split tool into tool and gamma
plot.df$gamma.char <- as.character(gsub("aldex._", "", plot.df$tool))
plot.df$gamma.num <- as.numeric(gsub("aldex._", "", plot.df$tool))
plot.df$tool <- gsub("_.$","", plot.df$tool)

plot.immuno <- plot.df %>% 
  filter(dataset == "immuno")


#################### Density Plot Panel ######################

group_means <- plot.immuno %>%
  group_by(tool, gamma.char) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

group_means0 <- plot.df %>%
  group_by(tool, dataset, gamma.num) %>%
  summarise(mean_value = mean(abs), .groups='drop') #####

p1<-plot.immuno %>% 
  ggplot(aes(x = abs, fill=gamma.char)) + 
  geom_density(alpha=0.5)+
  geom_vline(data = group_means, aes(xintercept = mean_value, colour=gamma.char), linetype="dashed", linewidth=0.5) +
  theme_bw()+
  facet_wrap(~tool)+
  xlim(0,2)+
  labs(title = "Immunotherapy Minimum Difference Density Distribution")+
  xlab('Minimum Absolute Difference for Significance')+
  ylab('Density')


p2 <- ggplot(group_means0, aes(x = gamma.num, y = mean_value, colour=dataset, shape = tool)) +
  geom_jitter(width = 0.01, alpha = 0.5) + 
  labs(title = "Average Minimum Difference Across Datasets",
       x = "Scale",
       y = "Mean Absolute Minimum Difference for Significance") +
  theme_bw()+
  geom_point(size = 2.5)

p1|p2 #to plot as a 2 panel figure

################# volcano plots ######################

load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

# Convert results to data frames
df1 <- as.data.frame(immuno.data_0.aldex2$t.data[[30]])
df2 <- as.data.frame(immuno.data_1.aldex2$t.data[[30]])
df3 <- as.data.frame(immuno.data_2.aldex2$t.data[[30]])
df4 <- as.data.frame(immuno.data_5.aldex2$t.data[[30]])

df1$scale <- "0"
df2$scale <- "0.1"  
df3$scale <- "0.2"
df4$scale <- "0.5" 

df_combined <- bind_rows(df1, df2, df3, df4)

# Volcano plot overlay
p3<- ggplot(df_combined, aes(x = diff.btw, y = -log10(we.eBH), color = scale)) +
       geom_point(alpha = 0.6) +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # FDR threshold
       theme_bw() +
       labs(title = "ALDEx2 Volcano Plot",
            x = "diff.btw",
            y = "-log10 Adjusted we.eBH",
            color = "scale")

# Convert results to data frames
df5 <- as.data.frame(immuno.data_0.aldex3$t.data[[30]])
df6 <- as.data.frame(immuno.data_1.aldex3$t.data[[30]])
df7 <- as.data.frame(immuno.data_2.aldex3$t.data[[30]])
df8 <- as.data.frame(immuno.data_5.aldex3$t.data[[30]])

df5$scale <- "0"
df6$scale <- "0.1"  
df7$scale <- "0.2"
df8$scale <- "0.5" 

df_combined2 <- bind_rows(df5, df6, df7, df8)

# Volcano plot overlay
p4<- ggplot(df_combined2, aes(x = estimate, y = -log10(p.val.adj), color = scale)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # FDR threshold
  theme_bw() +
  labs(title = "ALDEx3 Volcano Plot",
       x = "estimate",
       y = "-log10 Adjusted p.val.adj",
       color = "scale")


p3|p4 #to plot as 2 panel figure

#################### effect size plots #######################

#aldex2 panel
p5<-ggplot(df_combined, aes(x = diff.win, y = diff.btw, color = scale)) +
geom_point(alpha = 0.6, size = 0.8) +
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
theme_bw() +
labs(title = "ALDEx2 Effect Plot",
     x = "diff.win",
     y = "diff.btw",
     color = "scale")

p6<-ggplot(df_combined2, aes(x = std.error*sqrt(109), y = estimate, color = scale)) + #need to change this based on what corresponds to diff.win and diff.btw for aldex3
  geom_point(alpha = 0.6, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(title = "ALDEx3 Effect Plot",
       x = "std.error*sqrt(number of samples)",
       y = "estimate",
       color = "scale")

p5|p6



#plot(immuno.data_0.aldex3$t.data[[1]]$std.error*sqrt(120), immuno.data_0.aldex3$t.data[[1]]$estimate)
#std.error * sqrt(total number of samples) <- diff.win
#estimate <-  diff.btw

