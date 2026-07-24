# volcano & effect plots showing consequence of gamma: immunotherapy dataset

# Greg Gloor & Scott Dos Santos
# Last edited: 2026-07-24

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
              c(0,5), "aldex", rep(c("2.Rda", "3.Rda"), each = 2))

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
# repo.data <- paste0("~/Documents/GitHub/usri/data/imm_",
#                     paste0("gamma", rep(0:5, each = 2), "aldex", rep(c("2.Rda", "3.Rda"))))
# 
# 
# for(i in 1:length(ls(pattern = "gamma"))){
#   save(list = ls(pattern = "gamma")[i],
#        file = repo.data[i])
# }

########################### ALDEx3: data processing ###########################

# want to make volcano and effect plots for two different scenarios:
#
# 1) the old dual cutoff method (P <0.05 AND minimum 2-fold change) for genes 
# from thinned immuno data analysed when gamma = 0, but coloured based on 
# whether genes are significant at gamma = 0 AND when gamma = 0.5
#
# 2) using gamma only (i.e. no dual cutoff) for genes from thinned immuno data
# analysed when gamma = 0, but coloured based on whether genes are significant 
# at gamma = 0 AND when gamma = 0.5

# for scenario 1, will have grey lines showing -log10 P <0.05 and L2FC >1

# for scenario 2, will have grey lines showing -log10 P <0.05 and boundaries of
# significance at gamma = 0 and gamma = 0.5


# start by calculating -log10 BH-adjusted P values for both gamma value analyses
gamma0aldex3$log10p <- -log10(gamma0aldex3$p.val.adj)
gamma5aldex3$log10p <- -log10(gamma5aldex3$p.val.adj)

# calculate total number of Inf features and change to max value rounded to 
# nearest 2.5
length(which(gamma0aldex3$log10p == Inf)) # 155
max(gamma0aldex3$log10p[gamma0aldex3$log10p != Inf]) # 63.7
gamma0aldex3$log10p[which(gamma0aldex3$log10p == Inf)] <- 65

length(which(gamma5aldex3$log10p == Inf)) # 67
max(gamma5aldex3$log10p[gamma5aldex3$log10p != Inf]) # 44.7
gamma5aldex3$log10p[which(gamma5aldex3$log10p == Inf)] <- 45

# pull significant features for gamma = 0 and gamma = 0.5
a3.sig.g0 <- gamma0aldex3$p.val.adj <0.05
a3.sig.g5 <- gamma5aldex3$p.val.adj <0.05

# pull features which are ONLY significant when gamma = 0
a3.sig.g0 <- setdiff(which(a3.sig.g0), which(a3.sig.g5))

# pull features from gamma = 0 which are significant AND have L2FC >1
a3.sig.dual.more <- gamma0aldex3$p.val.adj <0.05 & abs(gamma0aldex3$estimate) >=1
a3.sig.dual.less <- gamma0aldex3$p.val.adj <0.05 & abs(gamma0aldex3$estimate) <1

# pull features which are not significant
a3.ns <- gamma0aldex3$p.val.adj >0.05

# subset aldex3 data at gamma = 0 for these five sets of features
plot.a3.sig.g0 <- gamma0aldex3[a3.sig.g0,]
plot.a3.sig.g5 <- gamma0aldex3[a3.sig.g5,]
plot.a3.sig.dualm <- gamma0aldex3[a3.sig.dual.more,]
plot.a3.sig.duall <- gamma0aldex3[a3.sig.dual.less,]
plot.a3.ns <- gamma0aldex3[a3.ns,]

# check all combinations add up to total non-Inf feature count
nrow(plot.a3.ns)+nrow(plot.a3.sig.duall)+nrow(plot.a3.sig.dualm) # 20,478
nrow(plot.a3.ns)+nrow(plot.a3.sig.g0)+nrow(plot.a3.sig.g5) # 20,478

# add title to aldex3 non-significant features df (edit this in inkscape)
plot.a3.ns$strip <- "ALDEx3: Plot title goes here"

########################## ALDEx3: dual cutoff plots ##########################

# volcano plot
cutoff.volc.a3 <- ggplot(data = plot.a3.ns, aes(x = estimate, y = log10p))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a3.sig.duall, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a3.sig.dualm, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = -1, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = 1, lty = 2, colour = "grey40", lwd = 0.5)+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("-Log"[10]~" P"["adj"]))+
  facet_wrap(~strip)+
  theme_bw()

cutoff.volc.a3


# effect plot
cutoff.eff.a3 <- ggplot(data = plot.a3.ns, aes(x = estimate, y = std.error))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a3.sig.duall, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a3.sig.dualm, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_vline(xintercept = -1, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = 1, lty = 2, colour = "grey40", lwd = 0.5)+
  scale_y_continuous(limits = c(0.04, 1.35), expand = c(0,0))+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("Standard error: Log"[2]~"fold change"))+
  facet_wrap(~strip)+
  theme_bw()

cutoff.eff.a3


# plot both together side by side
# png(paste0(repo.figs, "fig_a3_volcEffect_cutoff.png"),
#     units = "in", height = 4, width = 10, res = 600)

cutoff.volc.a3 | cutoff.eff.a3

# dev.off()

########################### ALDEx3: gamma only plots ###########################

# volcano plot
gamma.volc.a3 <- ggplot(data = plot.a3.ns, aes(x = estimate, y = log10p))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a3.sig.g0, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a3.sig.g5, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = -1, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = 1, lty = 2, colour = "grey40", lwd = 0.5)+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("-Log"[10]~" P"["adj"]))+
  facet_wrap(~strip)+
  theme_bw()

gamma.volc.a3


# effect plot
gamma.eff.a3 <- ggplot(data = plot.a3.ns, aes(x = estimate, y = std.error))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a3.sig.g0, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a3.sig.g5, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_abline(intercept = 0, slope = -0.3, lty = 3, colour = "grey40", lwd = 0.5)+
  geom_abline(intercept = 0, slope = 0.3, lty = 3, colour = "grey40", lwd = 0.5)+
  geom_abline(intercept = -0.21, slope = -0.3, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_abline(intercept = -0.21, slope = 0.3, lty = 2, colour = "grey40", lwd = 0.5)+
  scale_y_continuous(limits = c(0.04, 1.35), expand = c(0,0))+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("Standard error: Log"[2]~"fold change"))+
  facet_wrap(~strip)+
  theme_bw()

gamma.eff.a3


# plot both together side by side
# png(paste0(repo.figs, "fig_a3_volcEffect_gamma.png"),
#     units = "in", height = 4, width = 10, res = 600)

gamma.volc.a3 | gamma.eff.a3

# dev.off()

########################### ALDEx2: data processing ###########################

# want to make analagous volcano and effect plots for ALDEx2 because reviewers
# will probably ask for them

# start by calculating -log10 BH-adjusted P values for both gamma value analyses
gamma0aldex2$log10p <- -log10(gamma0aldex2$we.eBH)
gamma5aldex2$log10p <- -log10(gamma5aldex2$we.eBH)

# calculate total number of Inf features and change to max value rounded to 
# nearest 2.5
length(which(gamma0aldex2$log10p == Inf)) # 152
max(gamma0aldex2$log10p[gamma0aldex2$log10p != Inf]) # 58.9
gamma0aldex2$log10p[which(gamma0aldex2$log10p == Inf)] <- 60

length(which(gamma5aldex2$log10p == Inf)) # 50
max(gamma5aldex2$log10p[gamma5aldex2$log10p != Inf]) # 35.1
gamma5aldex2$log10p[which(gamma5aldex2$log10p == Inf)] <- 37.5

# pull significant features for gamma = 0 and gamma = 0.5
a2.sig.g0 <- gamma0aldex2$we.eBH <0.05
a2.sig.g5 <- gamma5aldex2$we.eBH <0.05

# pull features which are ONLY significant when gamma = 0
a2.sig.g0 <- setdiff(which(a2.sig.g0), which(a2.sig.g5))

# pull features from gamma = 0 which are significant AND have L2FC >1
a2.sig.dual.more <- gamma0aldex2$we.eBH <0.05 & abs(gamma0aldex2$diff.btw) >=1
a2.sig.dual.less <- gamma0aldex2$we.eBH <0.05 & abs(gamma0aldex2$diff.btw) <1

# pull features which are not significant
a2.ns <- gamma0aldex2$we.eBH >0.05

# subset aldex2 data at gamma = 0 for these five sets of features
plot.a2.sig.g0 <- gamma0aldex2[a2.sig.g0,]
plot.a2.sig.g5 <- gamma0aldex2[a2.sig.g5,]
plot.a2.sig.dualm <- gamma0aldex2[a2.sig.dual.more,]
plot.a2.sig.duall <- gamma0aldex2[a2.sig.dual.less,]
plot.a2.ns <- gamma0aldex2[a2.ns,]

# check all combinations add up to total non-Inf feature count
nrow(plot.a2.ns)+nrow(plot.a2.sig.duall)+nrow(plot.a2.sig.dualm) # 20,478
nrow(plot.a2.ns)+nrow(plot.a2.sig.g0)+nrow(plot.a2.sig.g5) # 20,478

# add title to aldex2 non-significant features df (edit this in inkscape)
plot.a2.ns$strip <- "ALDEx2: Plot title goes here"

########################## ALDEx2: dual cutoff plots ##########################

# volcano plot
cutoff.volc.a2 <- ggplot(data = plot.a2.ns, aes(x = diff.btw, y = log10p))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a2.sig.duall, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a2.sig.dualm, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = -1, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = 1, lty = 2, colour = "grey40", lwd = 0.5)+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("-Log"[10]~" P"["adj"]))+
  facet_wrap(~strip)+
  theme_bw()

cutoff.volc.a2


# effect plot
cutoff.eff.a2 <- ggplot(data = plot.a2.ns, aes(x = diff.btw, y = diff.win))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a2.sig.duall, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a2.sig.dualm, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_vline(xintercept = -1, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = 1, lty = 2, colour = "grey40", lwd = 0.5)+
  scale_y_continuous(limits = c(0.3, 12), expand = c(0,0))+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("Log"[2]~"dispersion"))+
  facet_wrap(~strip)+
  theme_bw()

cutoff.eff.a2


# plot both together side by side
# png(paste0(repo.figs, "fig_a2_volcEffect_cutoff.png"),
#     units = "in", height = 4, width = 10, res = 600)

cutoff.volc.a2 | cutoff.eff.a2

# dev.off()

########################### ALDEx2: gamma only plots ###########################

# volcano plot
gamma.volc.a2 <- ggplot(data = plot.a2.ns, aes(x = diff.btw, y = log10p))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a2.sig.g0, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a2.sig.g5, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = -1, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_vline(xintercept = 1, lty = 2, colour = "grey40", lwd = 0.5)+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("-Log"[10]~" P"["adj"]))+
  facet_wrap(~strip)+
  theme_bw()

gamma.volc.a2


# effect plot
gamma.eff.a2 <- ggplot(data = plot.a2.ns, aes(x = diff.btw, y = diff.win))+
  geom_point(colour = rgb(0,0,0,0.2), size = 1)+
  geom_point(data = plot.a2.sig.g0, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "red", alpha = 0.5)+
  geom_point(data = plot.a2.sig.g5, size = 1.5, stroke = 0.25,
             shape = 21, colour = "black", fill = "dodgerblue", alpha = 0.5)+
  geom_abline(intercept = 0.2, slope = -2, lty = 3, colour = "grey40", lwd = 0.5)+
  geom_abline(intercept = 0.2, slope = 2, lty = 3, colour = "grey40", lwd = 0.5)+
  geom_abline(intercept = -2.6, slope = -2.3, lty = 2, colour = "grey40", lwd = 0.5)+
  geom_abline(intercept = -3.9, slope = 2.8, lty = 2, colour = "grey40", lwd = 0.5)+
  scale_y_continuous(limits = c(0.3, 12), expand = c(0,0))+
  xlab(expression("Log"[2]~"fold change"))+
  ylab(expression("Log"[2]~"dispersion"))+
  facet_wrap(~strip)+
  theme_bw()

gamma.eff.a2


# plot both together side by side
# png(paste0(repo.figs, "fig_a2_volcEffect_gamma.png"),
#     units = "in", height = 4, width = 10, res = 600)

gamma.volc.a2 | gamma.eff.a2

# dev.off()
