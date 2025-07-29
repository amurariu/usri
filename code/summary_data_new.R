library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

source('code/summary.data.fun.R')

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
#edited plot starts here -----------------------------------------
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
