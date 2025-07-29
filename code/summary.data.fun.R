anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

sum.fun <- function(aldex2_0, aldex2_1, aldex2_2, aldex2_5, aldex3_0, aldex3_1, aldex3_2, aldex3_5){
#summary for gamma = 1e-3, ALDEx2
minmax2_0 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax2_0)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt0 <- aldex2_0$t.data[[i]]$diff.btw > 0
  lt0 <- aldex2_0$t.data[[i]]$diff.btw < 0
  sig <- aldex2_0$t.data[[i]]$we.eBH < 0.05
  minmax2_0[i,2] <- min(aldex2_0$t.data[[i]]$diff.btw[gt0 & sig])
  minmax2_0[i,1] <- max(aldex2_0$t.data[[i]]$diff.btw[lt0 & sig])
  absmin0[i] <- abs(min(aldex2_0$t.data[[i]]$diff.btw[gt0 & sig]))
  absmax0[i] <- abs(max(aldex2_0$t.data[[i]]$diff.btw[lt0 & sig]))
  } 
minmax2_0[,3] <- c(absmin0, absmax0)

  #summary for gamma = 0.1, ALDEx2
minmax2_1 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax2_1)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt1 <- aldex2_1$t.data[[i]]$diff.btw > 0
  lt1 <- aldex2_1$t.data[[i]]$diff.btw < 0
  sig1 <- aldex2_1$t.data[[i]]$we.eBH < 0.05
  minmax2_1[i,2] <- min(aldex2_1$t.data[[i]]$diff.btw[gt1 & sig1])
  minmax2_1[i,1] <- max(aldex2_1$t.data[[i]]$diff.btw[lt1 & sig1])
  absmin1[i] <- abs(min(aldex2_1$t.data[[i]]$diff.btw[gt1 & sig1]))
  absmax1[i] <- abs(max(aldex2_1$t.data[[i]]$diff.btw[lt1 & sig1]))
  } 
minmax2_1[,3] <- rbind(absmin1, absmax1)
  
  #summary for gamma = 0.2, ALDEx2
minmax2_2 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax2_2)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt2 <- aldex2_2$t.data[[i]]$diff.btw > 0
  lt2 <- aldex2_2$t.data[[i]]$diff.btw < 0
  sig2 <- aldex2_2$t.data[[i]]$we.eBH < 0.05
  minmax2_2[i,2] <- min(aldex2_2$t.data[[i]]$diff.btw[gt2 & sig2])
  minmax2_2[i,1] <- max(aldex2_2$t.data[[i]]$diff.btw[lt2 & sig2])
  absmin2[i] <- abs(min(aldex2_2$t.data[[i]]$diff.btw[gt2 & sig2]))
  absmax2[i] <- abs(max(aldex2_2$t.data[[i]]$diff.btw[lt2 & sig2]))
}
minmax2_2[,3] <- rbind(absmin2, absmax2)
  
  #summary for gamma = 0.5, ALDEx2
minmax2_5 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax2_5)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt5 <- aldex2_5$t.data[[i]]$diff.btw > 0
  lt5 <- aldex2_5$t.data[[i]]$diff.btw < 0
  sig5 <- aldex2_5$t.data[[i]]$we.eBH < 0.05
  minmax2_5[i,2] <- min(aldex2_5$t.data[[i]]$diff.btw[gt5 & sig5])
  minmax2_5[i,1] <- max(aldex2_5$t.data[[i]]$diff.btw[lt5 & sig5])
  absmin5[i] <- abs(min(aldex2_5$t.data[[i]]$diff.btw[gt5 & sig5]))
  absmax5[i] <- abs(max(aldex2_5$t.data[[i]]$diff.btw[lt5 & sig5]))
}
minmax2_5[,3] <- rbind(absmin5, absmax5)
  
  #summary for gamma = 0, ALDEx3
minmax3_0 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax3_0)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt3_0 <- aldex3_0$t.data[[i]]$estimate > 0
  lt3_0 <- aldex3_0$t.data[[i]]$estimate < 0
  sig3_0 <- aldex3_0$t.data[[i]]$p.val.adj < 0.05
  minmax3_0[i,2] <- min(aldex3_0$t.data[[i]]$estimate[gt3_0 & sig3_0])
  minmax3_0[i,1] <- max(aldex3_0$t.data[[i]]$estimate[lt3_0 & sig3_0])
  absmin3_0[i] <- abs(min(aldex3_0$t.data[[i]]$estimate[gt3_0 & sig3_0]))
  absmax3_0[i] <- abs(max(aldex3_0$t.data[[i]]$estimate[lt3_0 & sig3_0]))
}
minmax3_0[,3] <- rbind(absmin3_0, absmax3_0)
  
  #summary for gamma = 0.1, ALDEx3
minmax3_1 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax3_1)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt3_1 <- aldex3_1$t.data[[i]]$estimate > 0
  lt3_1 <- aldex3_1$t.data[[i]]$estimate < 0
  sig3_1 <- aldex3_1$t.data[[i]]$p.val.adj < 0.05
  minmax3_1[i,2] <- min(aldex3_1$t.data[[i]]$estimate[gt3_1 & sig3_1])
  minmax3_1[i,1] <- max(aldex3_1$t.data[[i]]$estimate[lt3_1 & sig3_1])
  absmin3_1[i] <- abs(min(aldex3_1$t.data[[i]]$estimate[gt3_1 & sig3_1]))
  absmax3_1[i] <- abs(max(aldex3_1$t.data[[i]]$estimate[lt3_1 & sig3_1]))
}
minmax3_1[,3] <- rbind(absmin3_1, absmax3_1)
  
  #summary for gamma = 0.2, ALDEx3
minmax3_2 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax3_2)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt3_2 <- aldex3_2$t.data[[i]]$estimate > 0
  lt3_2 <-  aldex3_2$t.data[[i]]$estimate < 0
  sig3_2 <-  aldex3_2$t.data[[i]]$p.val.adj < 0.05
  minmax3_2[i,2] <- min(aldex3_2$t.data[[i]]$estimate[gt3_2 & sig3_2])
  minmax3_2[i,1] <- max(aldex3_2$t.data[[i]]$estimate[lt3_2 & sig3_2])
  absmin3_2[i] <- abs(min(aldex3_2$t.data[[i]]$estimate[gt3_2 & sig3_2]))
  absmax3_2[i] <- abs(max(aldex3_2$t.data[[i]]$estimate[lt3_2 & sig3_2]))
}
minmax3_2[,3] <- rbind(absmin3_2, absmax3_2)
  
  #summary for gamma = 0.5, ALDEx3
minmax3_5 <- matrix(data=NA, ncol=3, nrow=200)
colnames(minmax3_5)<-c('max', 'min', 'abs')
for(i in 1:100){
  gt3_5 <- aldex3_5$t.data[[i]]$estimate > 0
  lt3_5 <- aldex3_5$t.data[[i]]$estimate < 0
  sig3_5 <- aldex3_5$t.data[[i]]$p.val.adj < 0.05
  minmax3_5[i,2] <- min(aldex3_5$t.data[[i]]$estimate[gt3_5 & sig3_5])
  minmax3_5[i,1] <- max(aldex3_5$t.data[[i]]$estimate[lt3_5 & sig3_5])
  absmin3_5[i] <- abs(min(aldex3_5$t.data[[i]]$estimate[gt3_5 & sig3_5]))
  absmax3_5[i] <- abs(max(aldex3_5$t.data[[i]]$estimate[lt3_5 & sig3_5]))
}
minmax3_5[,3] <- rbind(absmin3_5, absmax3_5)


return(list(minmax2_0 = minmax2_0, minmax2_1 = minmax2_1, minmax2_2 = minmax2_2, minmax2_5 = minmax2_5, minmax3_0 = minmax3_0, minmax3_1 = minmax3_1, minmax3_2 = minmax3_2, minmax3_5 = minmax3_5))
}