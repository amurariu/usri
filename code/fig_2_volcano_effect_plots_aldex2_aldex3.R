# representative volcano & effect plots using immuno dataset

# Andreea Murariu & Scott Dos Santos
# Last edited: 2025-08-14

library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# set path to directory containing all analysis objects
anal.path <- "~/Documents/GitHub/ext_analysis/"
repo.figs <- "~/Documents/GitHub/usri/figures/"

################################ volcano plots ################################

# load all aldex analysis objects from immuno dataset and extract an iteration
# (30th iteration, picked out of the air)
for(i in list.files(anal.path, pattern = "immuno.data.aldex")){
  load(paste0(anal.path,i))
  gamma <- gsub("aldex.*_","", str_split_i(i, "\\.", 3))
  tool <- gsub("_.*","", str_split_i(i, "\\.", 3))
  object <- paste0("immuno.data_", gamma, ".", tool)
  assign(x = paste0("gamma", gamma, tool),
         value = get(object)$t.data[[30]])
}

# remove temporary objects
rm(gamma,tool,object,i,
   list = ls(pattern = "immuno."))

# filter gamma = 0 for points where P <0.05 in gamma = 0 BUT the corresponding
# features are NOT significant in other gamma values (filtering like this only
# works because all data frames are in the same order and have the same number
# of rows)
a2.to.plot0 <- gamma0aldex2 %>% 
  filter(we.eBH <0.05,
         gamma1aldex2$we.eBH >=0.05,
         gamma2aldex2$we.eBH >=0.05,
         gamma3aldex2$we.eBH >=0.05,
         gamma4aldex2$we.eBH >=0.05,
         gamma5aldex2$we.eBH >=0.05)

# filter gamma = 0 separately for points where P <0.05 in gamma = 1 / 2 / 5
a2.to.plot1 <- gamma0aldex2 %>% 
  filter(gamma1aldex2$we.eBH <0.05)

a2.to.plot2 <- gamma0aldex2 %>% 
  filter(gamma2aldex2$we.eBH <0.05)

a2.to.plot3 <- gamma0aldex2 %>% 
  filter(gamma3aldex2$we.eBH <0.05)

a2.to.plot4 <- gamma0aldex2 %>% 
  filter(gamma4aldex2$we.eBH <0.05)

a2.to.plot5 <- gamma0aldex2 %>% 
  filter(gamma5aldex2$we.eBH <0.05)

# set colours for gamma values for ALDEx2 and ALDEx3 as in FDR/TPR figure
cols.volcano <- c("#E6E6FA","#CBE3FC","#B0E2FF","#66B3F6","#2171B5","#08306B",
                  brewer.pal(n = 9, name = "YlOrRd")[c(1,3,5,6,8,9)])

# cols.volcano <- c(brewer.pal(n = 9, name = "Blues")[c(3,5,7)], "navy",
#                   brewer.pal(n = 9, name = "OrRd")[c(3,5,7,9)])

# make legend colour vector for scale_colour_manual
leg.cols.a2 <- c("0" = cols.volcano[1], "0.1" = cols.volcano[2], "0.2" = cols.volcano[3],
                 "0.3" = cols.volcano[4], "0.4" = cols.volcano[5], "0.5" = cols.volcano[6])

# add column for strip title
gamma0aldex2$title <- "ALDEx2"

# plot volcano plots for all features in gamma = 0 and overlay points which
# are significant at increasing gamma values
# png(paste0(repo.figs, "fig2_volcanoPlotsALDEx2.png"),
#     units = "in", height = 5, width = 5, res = 600)

vol.a2 <- ggplot(data = gamma0aldex2, aes(x = diff.btw, y = -log10(we.eBH)))+
  geom_point(alpha = 0.4, size = 0.8, colour = "grey50")+
  geom_point(data = a2.to.plot0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+  # FDR threshold
  scale_y_continuous(limits = c(0,65), expand = c(0.005,0.75))+
  scale_x_continuous(limits = c(-7.75,6.25), expand = c(0.005, 0.005))+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a2)+
  labs(x = "Log difference between groups", y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 8), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "top")

vol.a2

# dev.off()


# same for ALDEx3: but first calculate dispersion (diff w/in groups), for the
# immuno dataset (which has 109 samples) analysed with ALDEx3
# diff.win = std.error * sqrt(no. samples)
gamma0aldex3 <- gamma0aldex3 %>% 
  mutate(diff.win = std.error * sqrt(109))

#filter gamma 0 
a3.to.plot0 <- gamma0aldex3 %>% 
  filter(p.val.adj <0.05,
         gamma1aldex3$p.val.adj >=0.05,
         gamma2aldex3$p.val.adj >=0.05,
         gamma3aldex3$p.val.adj >=0.05,
         gamma4aldex3$p.val.adj >=0.05,
         gamma5aldex3$p.val.adj >=0.05)

a3.to.plot1 <- gamma0aldex3 %>% 
  filter(gamma1aldex3$p.val.adj <0.05)

a3.to.plot2 <- gamma0aldex3 %>% 
  filter(gamma2aldex3$p.val.adj <0.05)

a3.to.plot3 <- gamma0aldex3 %>% 
  filter(gamma3aldex3$p.val.adj <0.05)

a3.to.plot4 <- gamma0aldex3 %>% 
  filter(gamma4aldex3$p.val.adj <0.05)

a3.to.plot5 <- gamma0aldex3 %>% 
  filter(gamma5aldex3$p.val.adj <0.05)

# legend colour vector
leg.cols.a3 <- c("0" = cols.volcano[7], "0.1" = cols.volcano[8], "0.2" = cols.volcano[9],
                 "0.3" = cols.volcano[10], "0.4" = cols.volcano[11], "0.5" = cols.volcano[12])

# strip title
gamma0aldex3$title <- "ALDEx3"

# plot
# png(paste0(repo.figs, "fig2_volcanoPlotsALDEx3.png"),
#     units = "in", height = 5, width = 5, res = 600)

vol.a3 <- ggplot(data = gamma0aldex3, aes(x = estimate, y = -log10(p.val.adj)))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = a3.to.plot0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+  
  geom_point(data = a3.to.plot5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+  # FDR threshold
  scale_y_continuous(limits = c(0,65), expand = c(0.005,0.75))+
  scale_x_continuous(limits = c(-7.75,6.25), expand = c(0.005, 0.005))+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a3)+
  labs(x = "Log difference between groups", y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 8), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "top")

vol.a3

# dev.off()


# both together on same plot
# png(paste0(repo.figs, "fig2_volcanoPlotsALDEx.png"),
#     units = "in", height = 5, width = 10, res = 600)

vol.a3.edit <- ggplot(data = gamma0aldex3, aes(x = estimate, y = -log10(p.val.adj)))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = a3.to.plot0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+  
  geom_point(data = a3.to.plot5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+  # FDR threshold
  scale_y_continuous(limits = c(0,65), expand = c(0.005,0.75))+
  scale_x_continuous(limits = c(-7.75,6.25), expand = c(0.005, 0.005))+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a3)+
  labs(x = "Log difference between groups", y = "")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 8), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "top",
        axis.title.y = element_blank())


vol.a2 | vol.a3.edit

# dev.off()

################################# effect plots #################################

# data for plotting difference within groups vs. difference between groups is
# already in all data frames: either natively, or calculated in previous section

# aldex2
# png(paste0(repo.figs, "fig2_effectPlotsALDEx2.png"),
#     units = "in", height = 5, width = 5, res = 600)

eff.a2 <- ggplot(data = gamma0aldex2, aes(x = diff.win, y = diff.btw))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = a2.to.plot0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a2.to.plot5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_abline(slope = 1,  intercept = 0, linetype = "dashed", color = "black")+
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black")+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a2)+
  scale_x_continuous(limits = c(0,14), expand = c(0.0001, 0.001))+
  scale_y_continuous(limits = c(-7.75,6.25), expand = c(0.001, 0.001))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 8), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "top")

eff.a2

# dev.off()

# aldex3
# png(paste0(repo.figs, "fig2_effectPlotsALDEx3.png"),
#     units = "in", height = 5, width = 5, res = 600)

eff.a3 <- ggplot(data = gamma0aldex3, aes(x = diff.win, y = estimate))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = a3.to.plot0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_abline(slope = 1,  intercept = 0, linetype = "dashed", color = "black")+
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black")+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a3)+
  scale_x_continuous(limits = c(0,14), expand = c(0.0001, 0.001))+
  scale_y_continuous(limits = c(-7.75,6.25), expand = c(0.001, 0.001))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 8), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "top")

eff.a3

# dev.off()

# both
# png(paste0(repo.figs, "fig2_effectPlotsALDEx.png"),
#     units = "in", height = 5, width = 10, res = 600)

eff.a3.edit <- ggplot(data = gamma0aldex3, aes(x = diff.win, y = estimate))+
  geom_point(alpha = 0.4, size = 0.9, colour = "grey50")+
  geom_point(data = a3.to.plot0, size = 1.5, aes(fill = "0"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot1, size = 1.5, aes(fill = "0.1"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot2, size = 1.5, aes(fill = "0.2"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot3, size = 1.5, aes(fill = "0.3"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot4, size = 1.5, aes(fill = "0.4"), shape = 21, colour = "black", stroke = 0.25)+
  geom_point(data = a3.to.plot5, size = 1.5, aes(fill = "0.5"), shape = 21, colour = "black", stroke = 0.25)+
  geom_abline(slope = 1,  intercept = 0, linetype = "dashed", color = "black")+
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black")+
  scale_fill_manual(name = "Scale (\u03b3)", values = leg.cols.a3)+
  scale_x_continuous(limits = c(0,14), expand = c(0.0001, 0.001))+
  scale_y_continuous(limits = c(-7.75,6.25), expand = c(0.001, 0.001))+
  labs(x = "Log difference within groups")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 8), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 9, face = "bold"), legend.position = "top",
        axis.title.y = element_blank())

eff.a2 | eff.a3.edit

# dev.off()
