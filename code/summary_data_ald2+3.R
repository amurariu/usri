anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

#for ALDEx2 with gamma = 1e-3
minmax2_0 <- matrix(data=NA, ncol=3, nrow=200)
for(i in 1:100){
  #summary for gamma = 0, ALDEx2
  gt0 <- immuno.data_0.aldex2$t.data[[i]]$diff.btw > 0
  lt0 <- immuno.data_0.aldex2$t.data[[i]]$diff.btw < 0
  sig <- immuno.data_0.aldex2$t.data[[i]]$we.eBH < 0.05
  minmax[i,2] <- min(immuno.data_0.aldex2$t.data[[i]]$diff.btw[gt0 & sig])
  minmax[i,1] <- max(immuno.data_0.aldex2$t.data[[i]]$diff.btw[lt0 & sig])
  absmin0 <- abs(min(immuno.data_0.aldex2$t.data[[i]]$diff.btw[gt0 & sig0]))
  absmax0 <- abs(max(immuno.data_0.aldex2$t.data[[i]]$diff.btw[lt0 & sig0]))
  }
  minmax[,3] <- rbind(absmin0, absmax0)
  
  minmax2_1 <- matrix(data=NA, ncol=3, nrow=200)
  for(i in 1:100){
    #summary for gamma = 0, ALDEx2
    gt1 <- immuno.data_1.aldex2$t.data[[i]]$diff.btw > 0
    lt1 <- immuno.data_1.aldex2$t.data[[i]]$diff.btw < 0
    sig1 <- immuno.data_1.aldex2$t.data[[i]]$we.eBH < 0.05
    minmax[i,2] <- min(immuno.data_1.aldex2$t.data[[i]]$diff.btw[gt0 & sig])
    minmax[i,1] <- max(immuno.data_1.aldex2$t.data[[i]]$diff.btw[lt0 & sig])
    absmin0 <- abs(min(immuno.data_1.aldex2$t.data[[i]]$diff.btw[gt0 & sig0]))
    absmax0 <- abs(max(immuno.data_1.aldex2$t.data[[i]]$diff.btw[lt0 & sig0]))
  }
  minmax[,3] <- rbind(absmin0, absmax0)