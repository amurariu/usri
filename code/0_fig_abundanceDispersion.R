# log abundance vs. dispersion for example dataset: ALDEx2 all datasets

# Scott Dos Santos
# Last edited: 2025-10-24

library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

# set path to GitHub repository
repo <- "~/Documents/GitHub/usri/"

################################# extract data #################################

# load in ALDEx2 output (gamma = 0) aggregated by dataset directly from GitHub
url.data <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/all_dispersionAbundance.Rda"
tf.data <- tempfile(fileext = ".Rda")
download.file(url.data, destfile = tf.data, mode = "wb")
load(tf.data)
unlink(tf.data)



# # original code for generating 'plot.df' loaded above
# 
# # set external analysis directory
# data <- "~/Documents/GitHub/ext_analysis/"
# 
# # pick random iteration of data to plot (reproducible)
# set.seed(1)
# inst <- sample(1:100, 1)
# 
# # make vector of datasets
# datasets <- c("brca","immuno","kirc","lihc","luad","mts","oral","prad","sccyto","thca")
# 
# # load aldex2 analysis objects from above datasets only for gamma = 0
# # and extract randomly selected instance (68 if set.seed = 1)
# tmplist <- list()
# for(i in datasets){
#   load(paste0(data,i, ".data.aldex2_0.Rda"))
#   tmpdf <- get(paste0(i, ".data_0.aldex2"))$t.data[[inst]]
#   tmpdf$dataset <- i
#   tmplist[[i]] <-tmpdf
#   rm(tmpdf)
# }
# 
# # add yeast data
# load(paste0(repo, "data/yeastUnf.A2g0.Rda"))
# yeastUnf.data_0.aldex2$dataset <- "yeast"
# tmplist[["yeast"]] <- yeastUnf.data_0.aldex2
# 
# rm(list = ls(pattern = "\\.aldex2$"))
# 
# # bind all dataframes
# plot.df <- do.call(rbind, tmplist)
# plot.df <- plot.df %>% 
#   relocate(dataset, .before = rab.all)
# 
# # save as compressed .Rda file
# # save(plot.df, file = paste0(repo, "data/all_dispersionAbundance.Rda"))
# 
# rm(tmplist)

######################## plot dispersion vs. abundance ########################

# # plot log relative abundance vs. dispersion for all datasets (may have to
# # reduce down to one bulk RNAseq, mts, yeast, scRNAseq)
# ggplot(data = plot.df, aes(x = diff.win, y = rab.all))+
#   geom_point(aes(colour = dataset), shape = 19, size = 0.75, alpha = 0.4)+
#   labs(x = "Log relative abundance", y = "Log dispersion within groups")+
#   theme_bw()

# needs filtering, based on above being pretty much illegible
plot.filt <- plot.df %>% 
  filter(dataset %in% c("immuno", "mts", "oral", "sccyto", "yeast")) %>% 
  mutate(strip = "ALDEx2: dispersion vs. relative abundance")

# make colour vector for datasets (n = 5)
cols.dataset <- c("purple2", "dodgerblue2", "olivedrab", "orangered2", "goldenrod2")

# make 4 separate graphs: main scatterplot in bottom left, density plots for x 
# and y values in top right and bottom right, then blank plot for top right
bl <- ggplot(data = plot.filt, aes(x = diff.win, y = rab.all))+
  geom_point(aes(colour = dataset), shape = 19, size = 0.2, alpha = 0.35, show.legend = F)+
  labs(x = "Log dispersion within groups", y = "Log relative abundance")+
  scale_colour_manual(values = cols.dataset)+
  theme_bw()+
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))

tl <- ggplot(data = plot.filt, aes(x = diff.win))+
  geom_density(aes(colour = dataset), linewidth = 0.4, show.legend = F)+
  scale_colour_manual(values = cols.dataset)+
  labs(y = "Density")+
  theme_bw()+
  facet_wrap(~strip)+
  theme(axis.title.x = element_blank(), plot.margin = margin(0.2,0.2,0,0.415,"cm"),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        axis.text.x = element_blank(), strip.text = element_text(face = "bold"))

br <- ggplot(data = plot.filt, aes(y = rab.all))+
  geom_density(aes(colour = dataset), linewidth = 0.4, show.legend = F)+
  scale_colour_manual(values = cols.dataset)+
  labs(x = "Density")+
  theme_bw()+
  theme(axis.title.y = element_blank(), plot.margin = margin(0.2,0.2,0.2,0,"cm"),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        axis.text.y = element_blank())

tr <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.x = element_blank(), 
        axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.line = element_blank())

# png(paste0(repo,"figures/supplFig_dispersionAbundance.png"),
#     units = "in", height = 4, width = 6, res = 600)

grid.arrange(tl,tr,bl,br, ncol = 2, nrow = 2, heights = c(1, 2), widths = c(4.25, 1))

# dev.off()


# make figure showing disp vs abund containing legend to use in inkscape later
# png(paste0(repo, "figures/supplFig_dispersionAbundance_legend.png"),
#     units = "in", height = 4, width = 6, res = 600)

bl.leg <- ggplot(data = plot.filt, aes(x = diff.win, y = rab.all))+
  geom_point(aes(colour = dataset), shape = 19, size = 0.75, alpha = 0.35)+
  labs(x = "Log dispersion within groups", y = "Log relative abundance")+
  scale_colour_manual(values = cols.dataset)+
  theme_bw()+
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))

bl.leg

# dev.off()
