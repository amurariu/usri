# representative volcano & effect plots using immunotherapy dataset

# Andreea Murariu & Scott Dos Santos
# Last edited: 2025-10-24

library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# path to figure directory on GH
repo.figs <- "~/Documents/GitHub/usri/figures/"

#################################### setup ####################################

# load in all ALDEx2 and ALDEx3 results at gamma = 0 to gamma = 0.5 for the 
# immunotherapy dataset's 30th iteration (this number was chosen out of the air)

# files to load
rda <- paste0("https://github.com/amurariu/usri/raw/refs/heads/main/data/imm_gamma",
              0:5, "aldex", rep(c("2.Rda", "3.Rda"), each = 6))

# load directly from GitHub
for(i in rda){
  tf <- tempfile(fileext = ".Rda")
  download.file(url = i, destfile = tf, mode = "wb")
  load(tf)
  unlink(tf)
  rm(tf)
}



# # original code to produce all 'gammaXaldexY' data frames
# 
# # load all aldex analysis objects from immuno dataset and extract the 30th
# # iteration- picked out of the air (these loops are horrible to read- sorry to
# # whoever lays eyes on them)
# 
# # path to directory containing all analysis objects
# anal.path <- "~/Documents/GitHub/ext_analysis/"
# 
# # load analysis outputs and extract 30th loops if they don't already exist
# # NOTE: requires generation of analysis objects to run from scratch!
# for(i in list.files(anal.path, pattern = "immuno.data.aldex")){
#   load(paste0(anal.path,i))
# 
#   gamma <- gsub("aldex.*_","", str_split_i(i, "\\.", 3))
#   tool <- gsub("_.*","", str_split_i(i, "\\.", 3))
#   object <- paste0("immuno.data_", gamma, ".", tool)
# 
#   assign(x = paste0("gamma", gamma, tool),
#          value = get(object)$t.data[[30]])
# }
# 
# # remove temporary objects
# rm(gamma,tool,object,i, list = ls(pattern = "immuno."))
# 
# # write 30th iteration of immuno data to files
# repo.data <- paste("~/Documents/GitHub/usri/data/imm_",
#                    paste0("gamma", 0:5, "aldex", rep(c("2.Rda", "3.Rda"), each = 6)), sep = "")
# 
# 
# for(i in 1:length(ls(pattern = "gamma"))){
#   save(list = ls(pattern = "gamma")[i],
#        file = repo.data[i])
# }

################################ volcano plots ################################

# calculate -log10 P values, then change any Inf to reasonable value to be
# displayed on plot (62.5 for A2)
gamma0aldex2$qval <- -log10(gamma0aldex2$we.eBH)
gamma0aldex2$qval <- case_when(gamma0aldex2$qval == Inf ~ 62.5,
                               .default = gamma0aldex2$qval)

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
# png(paste0(repo.figs, "fig_volcanoPlotsALDEx2.png"),
#     units = "in", height = 5, width = 5, res = 600)

vol.a2 <- ggplot(data = gamma0aldex2, aes(x = diff.btw, y = qval))+
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
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10))

vol.a2

# dev.off()


# same for ALDEx3: but first calculate dispersion (diff w/in groups), for the
# immuno dataset (which has 109 samples) analysed with ALDEx3
# diff.win = std.error * sqrt(no. samples)
gamma0aldex3 <- gamma0aldex3 %>% 
  mutate(diff.win = std.error * sqrt(109),
         qval = -log10(p.val.adj))

gamma0aldex3$qval <- case_when(gamma0aldex3$qval == Inf ~ 62.5,
                               .default = gamma0aldex3$qval)

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
# png(paste0(repo.figs, "fig_volcanoPlotsALDEx3.png"),
#     units = "in", height = 5, width = 5, res = 600)

vol.a3 <- ggplot(data = gamma0aldex3, aes(x = estimate, y = qval))+
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
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10))

vol.a3

# dev.off()


# both together on same plot
# png(paste0(repo.figs, "fig_volcanoPlotsALDEx.png"),
#     units = "in", height = 5, width = 10, res = 600)

vol.a3.edit <- ggplot(data = gamma0aldex3, aes(x = estimate, y = qval))+
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
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10),
        axis.title.y = element_blank())

vol.a2 | vol.a3.edit

# dev.off()

################################# effect plots #################################

# data for plotting difference within groups vs. difference between groups is
# already in all data frames: either natively, or calculated in previous section

# aldex2
# png(paste0(repo.figs, "fig_effectPlotsALDEx2.png"),
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
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10))

eff.a2

# dev.off()

# aldex3
# png(paste0(repo.figs, "fig_effectPlotsALDEx3.png"),
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
  labs(x = "Log error of linear model", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10))

eff.a3

# dev.off()

# both
# png(paste0(repo.figs, "fig_effectPlotsALDEx.png"),
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
  labs(x = "Log error of linear model")+
  theme_bw()+
  facet_wrap(~title)+
  guides(fill = guide_legend(nrow = 1))+
  theme(legend.box.spacing = unit(0.01, "cm"), legend.key.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = 9), strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 10, face = "bold"), legend.position = "top",
        axis.title = element_text(size = 10), axis.text = element_text(size = 10),
        axis.title.y = element_blank())

eff.a2 | eff.a3.edit

# dev.off()

########################### calculating sig features ###########################

# load in number of significantly different genes at each ALDEx2 or ALDEx3
# iteration for the immunotherapy dataset at each gamma value
url.sig <- "https://github.com/amurariu/usri/raw/refs/heads/main/data/imm_sigGenesALDEx.Rda"
tf.sig <- tempfile(fileext = ".Rda")
download.file(url = url.sig, destfile = tf.sig, mode = "wb")
load(tf.sig)
unlink(tf.sig)
rm(list = c("url.sig", "tf.sig"))

# make df and vectors for holding data
sig.genes.sum <- data.frame(matrix(data = NA, nrow = 6, ncol = 4, 
                                   dimnames = list(c(paste0("gamma", seq(0,5,1))),
                                                   c("A2mean","A2sd", "A3mean", "A3sd"))))

# calculate means and standard deviations for ALDEx2 and ALDEx3
sig.genes.sum[,1] <- apply(sig.genes.all[,1:6], 2, mean)
sig.genes.sum[,2] <- apply(sig.genes.all[,1:6], 2, sd)
sig.genes.sum[,3] <- apply(sig.genes.all[,7:12], 2, mean)
sig.genes.sum[,4] <- apply(sig.genes.all[,7:12], 2, sd)



# # original code for producing 'sig.genes'
# 
# # make vectors for holding number of significant genes at each iteration
# for(i in 0:5){
#   assign(x = paste0("sig", i, ".a2"), value = vector())
#   assign(x = paste0("sig", i, ".a3"), value = vector())
# }
# 
# # loop over all aldex outputs for immuno data and calculate the number of 
# # significantly different features at each gamma value for ALDEx2 and ALDEx3
# # separately
# for(i in 1:100){
#   
#   # pull ALDEx2 & ALDEx3 data for current iteration
#   for(j in 0:5){
#     assign(x = paste0("a2.df.", j),
#            value = get(paste0("immuno.data_", j, ".aldex2"))$t.data[[i]])
#     
#     assign(x = paste0("a3.df.", j),
#            value = get(paste0("immuno.data_", j, ".aldex3"))$t.data[[i]])
#     }
#   
#   # get the number of significant features for each scale value
#   for(j in 0:5){
#     assign(x = paste0("sig", j, ".a2"),
#            value = append(x = get(paste0("sig", j, ".a2")), after = i-1, 
#                           values = nrow(get(paste0("a2.df.", j)) %>% filter(we.eBH <0.05))))
#     
#     assign(x = paste0("sig", j, ".a3"),
#            value = append(x = get(paste0("sig", j, ".a3")), after = i-1, 
#                           values = nrow(get(paste0("a3.df.", j)) %>% filter(p.val.adj <0.05))))
#     
#     rm(list = c(paste0("a2.df.", j), paste0("a3.df.", j)))
#     
#     }
# }
# 
# # bind vectors and save object
# sig.genes.all <- cbind(sig0.a2, sig1.a2, sig2.a2, sig3.a2, sig4.a2, sig5.a2,
#                       sig0.a3, sig1.a3, sig2.a3, sig3.a3, sig4.a3, sig5.a3)
# 
# save(sig.genes.all, file = paste0("~/Documents/GitHub/usri/data/imm_sigGenesALDEx.Rda"))
#   
# # delete temp files
# for(i in 1:5){
#   rm(list = c(paste0("sig",i,".a2"), paste0("sig",i,".a3")))
# }
