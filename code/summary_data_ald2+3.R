anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

minmax <- matrix(data=NA, ncol=2, nrow=100)
for(i in 16:100){
  gt0 <- immuno.data_0.aldex2$t.data[[i]]$diff.btw > 0
  lt0 <- immuno.data_0.aldex2$t.data[[i]]$diff.btw < 0
  sig0 <- immuno.data_0.aldex2$t.data[[i]]$we.eBH < 0.05
  minmax[i,2] <- min(immuno.data_0.aldex2$t.data[[i]]$diff.btw[gt0 & sig0])
  minmax[i,1] <- max(immuno.data_0.aldex2$t.data[[i]]$diff.btw[lt0 & sig0])
  
  gt1 <- immuno.data_1.aldex2$t.data[[i]]$diff.btw > 0
  lt1 <- immuno.data_1.aldex2$t.data[[i]]$diff.btw < 0
  sig1 <- immuno.data_1.aldex2$t.data[[i]]$we.eBH < 0.05
  minmax[i,4] <- min(immuno.data_1.aldex2$t.data[[i]]$diff.btw[gt1 & sig1])
  minmax[i,3] <- max(immuno.data_1.aldex2$t.data[[i]]$diff.btw[lt1 & sig1])
  
  gt2 <- immuno.data_2.aldex2$t.data[[i]]$diff.btw > 0
  lt2 <- immuno.data_2.aldex2$t.data[[i]]$diff.btw < 0
  sig2 <- immuno.data_2.aldex2$t.data[[i]]$we.eBH < 0.05
  minmax[i,6] <- min(immuno.data_2.aldex2$t.data[[i]]$diff.btw[gt2 & sig2])
  minmax[i,5] <- max(immuno.data_2.aldex2$t.data[[i]]$diff.btw[lt2 & sig2])
  
  gt5 <- immuno.data_5.aldex2$t.data[[i]]$diff.btw > 0
  lt5 <- immuno.data_5.aldex2$t.data[[i]]$diff.btw < 0
  sig5 <- immuno.data_5.aldex2$t.data[[i]]$we.eBH < 0.05
  minmax[i,8] <- min(immuno.data_5.aldex2$t.data[[i]]$diff.btw[gt5 & sig5])
  minmax[i,7] <- max(immuno.data_5.aldex2$t.data[[i]]$diff.btw[lt5 & sig5])
  
  gt3_0 <- immuno.data_0.aldex3$t.data[[i]]$estimate > 0
  lt3_0 <- immuno.data_0.aldex3$t.data[[i]]$estimate < 0
  sig3_0 <- immuno.data_0.aldex3$t.data[[i]]$we.eBH < 0.05
  minmax[i,10] <- min(immuno.data_0.aldex3$t.data[[i]]$estimate[gt3_0 & sig3_0])
  minmax[i,9] <- max(immuno.data_0.aldex3$t.data[[i]]$estimate[lt3_0 & sig3_0])
  
  gt3_1 <- immuno.data_1.aldex3$t.data[[i]]$estimate > 0
  lt3_1 <- immuno.data_1.aldex3$t.data[[i]]$estimate < 0
  sig3_0 <- immuno.data_0.aldex3$t.data[[i]]$we.eBH < 0.05
  minmax[i,12] <- min(immuno.data_1.aldex3$t.data[[i]]$estimate[gt3_1 & sig3_1])
  minmax[i,11] <- max(immuno.data_1.aldex3$t.data[[i]]$estimate[lt3_1 & sig3_1])
  
  gt3_2 <- immuno.data_2.aldex3$t.data[[i]]$estimate > 0
  lt3_2 <- immuno.data_2.aldex3$t.data[[i]]$estimate < 0
  sig3_2 <- immuno.data_2.aldex3$t.data[[i]]$we.eBH < 0.05
  minmax[i,14] <- min(immuno.data_2.aldex3$t.data[[i]]$estimate[gt3_2 & sig3_2])
  minmax[i,13] <- max(immuno.data_2.aldex3$t.data[[i]]$estimate[lt3_2 & sig3_2])
  
  gt3_5 <- immuno.data_5.aldex3$t.data[[i]]$estimate > 0
  lt3_5 <- immuno.data_5.aldex3$t.data[[i]]$estimate < 0
  sig3_5 <- immuno.data_5.aldex3$t.data[[i]]$we.eBH < 0.05
  minmax[i,16] <- min(immuno.data_5.aldex3$t.data[[i]]$estimate[gt3_5 & sig3_5])
  minmax[i,15] <- max(immuno.data_5.aldex3$t.data[[i]]$estimate[lt3_5 & sig3_5])
}

