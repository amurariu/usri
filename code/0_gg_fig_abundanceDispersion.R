

load(url('https://github.com/amurariu/usri/blob/main/data/ss.bulk_x.all.Rda?raw=TRUE'))

# need to label rows as ssbulk.1:d
# then have first column "ssbulk"

# then rbind to data in 0_fig_abundanceDispersion.R and follow through with
# new dataset to plot

names <- paste("ssbulk", 1:nrow(ss.bulk_x.all), sep=".")
dataset <- rep("ssbulk", length(names))
ssbulk <- data.frame(dataset, ss.bulk_x.all) 
rownames(ssbulk) <- names
colnames(ssbulk) <- c("dataset","rab.all","rab.win.0","rab.win.1","diff.btw","diff.win","effect","overlap","we.ep","we.eBH","wi.ep","wi.eBH")

plot.df <- rbind(plot.df, ssbulk)

### modified existing
## I fucking give up, how do you re-order points in ggplot without
## yet another package in the tidyverse -- that changes every release
# Scott - I want the ssbulk to be plotted at the lowest layer so that the sccyto and tiyani are overplotted on it
plot.filt <- plot.df %>% 
  filter(dataset %in% c("immuno", "mts", "ssbulk", "tiyani-16s", "sccyto", "yeast")) %>% 
  mutate(strip = "ALDEx2: dispersion vs. relative abundance")

plot.filt$dataset <- gsub("ssbulk", "bulk", plot.filt$dataset)

# calculate median dispersion values per dataset for overlaying on density plot
grp.meds <- plot.filt %>% 
  group_by(dataset) %>% 
  summarise(med = median(diff.win))

# make colour vector for datasets (n = 5), corresponding to line graph of
# mean minimum log differences by gamma
cols.dataset <- c("cyan","darkolivegreen4", "dodgerblue2", "orangered2", "grey40", "goldenrod2")

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
  geom_vline(data = grp.meds, aes(xintercept = med, colour = dataset),
             linetype = 2, linewidth = 0.3, show.legend = F)+
  scale_colour_manual(values = cols.dataset)+
  scale_y_continuous(limits = c(0,4), expand = c(0,0))+
  labs(y = "Density")+
  theme_bw()+
  facet_wrap(~strip)+
  theme(axis.title.x = element_blank(), plot.margin = margin(0.2,0.2,0,0.415,"cm"),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        axis.text.x = element_blank(), strip.text = element_text(face = "bold"))

br <- ggplot(data = plot.filt, aes(y = rab.all))+
  geom_density(aes(colour = dataset), linewidth = 0.4, show.legend = F)+
  scale_colour_manual(values = cols.dataset)+
  scale_x_continuous(limits = c(0,4), expand = c(0,0))+
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
