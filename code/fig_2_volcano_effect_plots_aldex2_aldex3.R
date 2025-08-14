#ALDEx2: need to keep columns diff.btw, we.ebH, diff.win
#ALDEx3: need to keep columns estimate, p.adj.val, std. error
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(RColorBrewer)
library(ALDEx2)

anal.path <- "~/Documents/GitHub/ext_analysis/"

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
df1.win <- as.data.frame(immuno.data_0.aldex2$t.data[[30]]$diff.win)
colnames(df1.win) <- c("diff.win")
df1.padj <- as.data.frame(immuno.data_0.aldex2$t.data[[30]]$we.eBH)
colnames(df1.padj) <- c("p.adj")
df1.btw <- as.data.frame(immuno.data_0.aldex2$t.data[[30]]$diff.btw)
colnames(df1.btw) <- c("diff.btw")
df1<-bind_cols(df1.win, df1.padj, df1.btw)
df1$scale <- "0"

df2.win <- as.data.frame(immuno.data_1.aldex2$t.data[[30]]$diff.win)
colnames(df2.win) <- c("diff.win")
df2.padj <- as.data.frame(immuno.data_1.aldex2$t.data[[30]]$we.eBH)
colnames(df2.padj) <- c("p.adj")
df2.btw <- as.data.frame(immuno.data_1.aldex2$t.data[[30]]$diff.btw)
colnames(df2.btw) <- c("diff.btw")
df2<-bind_cols(df2.win, df2.padj, df2.btw)
df2$scale <- "0.1"  

df3.win <- as.data.frame(immuno.data_2.aldex2$t.data[[30]]$diff.win)
colnames(df3.win) <- c("diff.win")
df3.padj <- as.data.frame(immuno.data_2.aldex2$t.data[[30]]$we.eBH)
colnames(df3.padj) <- c("p.adj")
df3.btw <- as.data.frame(immuno.data_2.aldex2$t.data[[30]]$diff.btw)
colnames(df3.btw) <- c("diff.btw")
df3<-bind_cols(df3.win, df3.padj, df3.btw)
df3$scale <- "0.2"

df4.win <- as.data.frame(immuno.data_5.aldex2$t.data[[30]]$diff.win)
colnames(df4.win) <- c("diff.win")
df4.padj <- as.data.frame(immuno.data_5.aldex2$t.data[[30]]$we.eBH)
colnames(df4.padj) <- c("p.adj")
df4.btw <- as.data.frame(immuno.data_5.aldex2$t.data[[30]]$diff.btw)
colnames(df4.btw) <- c("diff.btw")
df4<-bind_cols(df4.win, df4.padj, df4.btw)
df4$scale <- "0.5" 

df_combined <- bind_rows(df1, df2, df3, df4)
df_combined$tool <- "ALDEx2"


#ALDEx3 new data frames
df5.win <- as.data.frame(immuno.data_0.aldex3$t.data[[30]]$std.error*sqrt(109))
colnames(df5.win) <- c("diff.win")
df5.padj <- as.data.frame(immuno.data_0.aldex3$t.data[[30]]$p.val.adj)
colnames(df5.padj) <- c("p.adj")
df5.btw <- as.data.frame(immuno.data_0.aldex3$t.data[[30]]$estimate)
colnames(df5.btw) <- c("diff.btw")
df5<-bind_cols(df5.win, df5.padj, df5.btw)
df5$scale <- "0"

df6.win <- as.data.frame(immuno.data_1.aldex3$t.data[[30]]$std.error*sqrt(109))
colnames(df6.win) <- c("diff.win")
df6.padj <- as.data.frame(immuno.data_1.aldex3$t.data[[30]]$p.val.adj)
colnames(df6.padj) <- c("p.adj")
df6.btw <- as.data.frame(immuno.data_1.aldex3$t.data[[30]]$estimate)
colnames(df6.btw) <- c("diff.btw")
df6<-bind_cols(df6.win, df6.padj, df6.btw)
df6$scale <- "0.1"  

df7.win <- as.data.frame(immuno.data_2.aldex3$t.data[[30]]$std.error*sqrt(109))
colnames(df7.win) <- c("diff.win")
df7.padj <- as.data.frame(immuno.data_2.aldex3$t.data[[30]]$p.val.adj)
colnames(df7.padj) <- c("p.adj")
df7.btw <- as.data.frame(immuno.data_2.aldex3$t.data[[30]]$estimate)
colnames(df7.btw) <- c("diff.btw")
df7<-bind_cols(df7.win, df7.padj, df7.btw)
df7$scale <- "0.2"

df8.win <- as.data.frame(immuno.data_5.aldex3$t.data[[30]]$std.error*sqrt(109))
colnames(df8.win) <- c("diff.win")
df8.padj <- as.data.frame(immuno.data_5.aldex3$t.data[[30]]$p.val.adj)
colnames(df8.padj) <- c("p.adj")
df8.btw <- as.data.frame(immuno.data_5.aldex3$t.data[[30]]$estimate)
colnames(df8.btw) <- c("diff.btw")
df8<-bind_cols(df8.win, df8.padj, df8.btw)
df8$scale <- "0.5" 

df2_combined <- bind_rows(df5, df6, df7, df8)
df2_combined$tool <- "ALDEx3"

total <- bind_rows(df_combined, df2_combined)


# Volcano plot
p3<- ggplot(total, aes(x = diff.btw, y = -log10(p.adj), color = scale)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # FDR threshold
  theme_bw() +
  facet_wrap(~tool)+
  labs(title = "ALDEx2 Volcano Plot",
       x = "diff.btw",
       y = "-log10 Adjusted p.adj",
       color = "scale")



############ volcano plot with other colours: #####################
library(ggplot2)
library(RColorBrewer)

# 4 shades of blue for ALDEx2, 4 shades of orange for ALDEx3
cols.ald2 <- brewer.pal(9, "Blues")[4:7]
cols.ald3 <- brewer.pal(9, "Oranges")[4:7]

# Named vector matching interaction(tool, scale)
cols.dataset <- c(
  setNames(cols.ald2, paste0("ALDEx2.", unique(total$scale))),
  setNames(cols.ald3, paste0("ALDEx3.", unique(total$scale)))
)

p3 <- ggplot(total, aes(x = diff.btw, y = -log10(p.adj), color = interaction(tool, scale))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(name = "Tool & Scale", values = cols.dataset) +
  theme_bw() +
  facet_wrap(~tool) +
  labs(title = "Volcano Plot", x = "diff.btw", y = "-log10 Adjusted p.adj")


#################### effect size plots #######################

cols.ald2 <- brewer.pal(n = 9, name = "Blues")[4:7]
cols.ald3 <- brewer.pal(n = 9, name = "Oranges")[4:7]

cols.dataset <- c(
  setNames(cols.ald2, paste0("ALDEx2.", unique(total$scale))),
  setNames(cols.ald3, paste0("ALDEx3.", unique(total$scale)))
)

ggplot(total, aes(x = diff.win, y = diff.btw)) +
  geom_point(data = subset(total, p.adj > 0.05), color = "gray80", alpha = 0.3, size = 1) +
  geom_point(data = subset(total, p.adj <= 0.05), aes(color = interaction(tool, scale)), alpha = 0.8, size = 1) +
  scale_colour_manual(name = "Tool & Scale", values = cols.dataset) +
  geom_abline(slope = 1,  intercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  facet_wrap(scale ~ tool, ncol = 2) +
  labs(title = "Effect-Size Plot", 
       x = "diff.win", 
       y = "diff.btw")


#plot from figure
par(mfrow=c(1,1))
aldex.plot(immuno.data_0.aldex2$t.data[[10]], type='MW')
title("D: Effect", adj=0, line= 0.8) 
points(immuno.data_0.aldex2$t.data[[10]]$diff.win[which(immuno.data_2.aldex2$t.data[[10]]$we.eBH < 0.05)], 
       (immuno.data_0.aldex2$t.data[[10]]$diff.btw[which(immuno.data_2.aldex2$t.data[[10]]$we.eBH < 0.05)]), cex=0.6, pch=19, 
       col='orange')
points(immuno.data_0.aldex2$t.data[[10]]$diff.win[which(immuno.data_5.aldex2$t.data[[10]]$we.eBH < 0.05)], 
       (immuno.data_0.aldex2$t.data[[10]]$diff.btw[which(immuno.data_5.aldex2$t.data[[10]]$we.eBH < 0.05)]), col='blue', 
       cex=0.6, pch=19)

abline(v=c(-0.5,0.5), lty=2, col='darkgrey')
abline(v=c(-1.5,1.5), lty=2, col='darkgrey')


#plot(immuno.data_0.aldex3$t.data[[1]]$std.error*sqrt(120), immuno.data_0.aldex3$t.data[[1]]$estimate)
#std.error * sqrt(total number of samples) <- diff.win
#estimate <-  diff.btw


#test code
cols.ald2 <- brewer.pal(n = 9, name = "Blues")[4:7]
cols.ald3 <- brewer.pal(n = 9, name = "Oranges")[4:7]

cols.dataset <- c(
  setNames(cols.ald2, paste0("ALDEx2.", unique(total$scale))),
  setNames(cols.ald3, paste0("ALDEx3.", unique(total$scale)))
)


gamma0 <- immuno.data_0.aldex2$t.data[[20]]
gamma2 <- immuno.data_2.aldex2$t.data[[20]]
gamma5 <- immuno.data_5.aldex2$t.data[[20]]

gamma0$title <- "ALDEx2"

to.plot.0 <- gamma0 %>%
  filter(we.eBH < 0.05, immuno.data_2.aldex2$t.data[[20]]$we.eBH >= 0.05,
         immuno.data_5.aldex2$t.data[[20]]$we.eBH >= 0.05)

to.plot.2 <- gamma2 %>%
  filter(we.eBH < 0.05,  immuno.data_5.aldex2$t.data[[20]]$we.eBH >= 0.05)

to.plot.5 <- gamma5 %>%
  filter(immuno.data_5.aldex2$t.data[[20]]$we.eBH < 0.05)

##newest version
ggplot(gamma0, aes(diff.win, diff.btw)) +
geom_point(alpha = 0.6, size = 1)+
geom_point(data = to.plot.0, color = "red", size = 1.5)+
geom_point(data = to.plot.2, color = "orange", size = 1.5)+
  geom_point(data = to.plot.5, color = "blue", size = 1.5)+
  geom_abline(slope = 1,  intercept = 0, linetype = "dashed", color = "black")+
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black")+
  labs(title = "D: Effect-Size Plots", x = "diff.win", y = "diff.btw")+
  theme_bw()



gamma3.0 <- immuno.data_0.aldex3$t.data[[20]]
gamma3.2 <- immuno.data_2.aldex3$t.data[[20]]
gamma3.5 <- immuno.data_5.aldex3$t.data[[20]]

gamma0$title <- "ALDEx2"

to.plot.3.0 <- gamma3.0 %>%
  filter(p.val.adj < 0.05, immuno.data_2.aldex3$t.data[[20]]$p.val.adj >= 0.05,
         immuno.data_5.aldex3$t.data[[20]]$p.val.adj >= 0.05)

to.plot.3.2 <- gamma3.2 %>%
  filter(p.val.adj < 0.05,  immuno.data_5.aldex3$t.data[[20]]$p.val.adj >= 0.05)

to.plot.3.5 <- gamma3.5 %>%
  filter(immuno.data_5.aldex3$t.data[[20]]$p.val.adj < 0.05)

##newest version
ggplot(gamma0, aes(diff.win, diff.btw)) +
  geom_point(alpha = 0.6, size = 1)+
  geom_point(data = to.plot.0, color = "red", size = 1.5)+
  geom_point(data = to.plot.2, color = "orange", size = 1.5)+
  geom_point(data = to.plot.5, color = "blue", size = 1.5)+
  geom_abline(slope = 1,  intercept = 0, linetype = "dashed", color = "black")+
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "black")+
  labs(title = "D: Effect-Size Plots", x = "diff.win", y = "diff.btw")+
  theme_bw()




