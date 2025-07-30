library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

anal.path <- "../ext_analysis/"
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

