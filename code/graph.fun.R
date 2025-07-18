source('code/get_confusion.R')
#DESeq_var, edgeR_var, limma_var,

graph.fun <- function(aldex_0_var, aldex_1_var, aldex_2_var, aldex_5_var, aldex3_0_var, aldex3_1_var, aldex3_2_var, aldex3_5_var){
  
  #perform get_confusion function first
  # des.conf <- get_confusion(DESeq_var, "DESeq")
  # edge.conf <- get_confusion(edgeR_var, "edgeR")
  ald0.conf <- get_confusion(aldex_0_var, "ALDEx2", FDR=0.1)
  ald1.conf <- get_confusion(aldex_1_var, "ALDEx2", FDR=0.1)
  ald2.conf <- get_confusion(aldex_2_var, "ALDEx2", FDR=0.1)
  ald5.conf <- get_confusion(aldex_5_var, "ALDEx2", FDR=0.1)
  ald3_0.conf <- get_confusion(aldex3_0_var, "ALDEx3", FDR=0.1)
  ald3_1.conf <- get_confusion(aldex3_1_var, "ALDEx3", FDR=0.1)
  ald3_2.conf <- get_confusion(aldex3_2_var, "ALDEx3", FDR=0.1)
  ald3_5.conf <- get_confusion(aldex3_5_var, "ALDEx3", FDR=0.1)
  ald0.conf5 <- get_confusion(aldex_0_var, "ALDEx2", FDR=0.05)
  ald1.conf5 <- get_confusion(aldex_1_var, "ALDEx2", FDR=0.05)
  ald2.conf5 <- get_confusion(aldex_2_var, "ALDEx2", FDR=0.05)
  ald5.conf5 <- get_confusion(aldex_5_var, "ALDEx2", FDR=0.05)
  ald3_0.conf5 <- get_confusion(aldex3_0_var, "ALDEx3", FDR=0.05)
  ald3_1.conf5 <- get_confusion(aldex3_1_var, "ALDEx3", FDR=0.05)
  ald3_2.conf5 <- get_confusion(aldex3_2_var, "ALDEx3", FDR=0.05)
  ald3_5.conf5 <- get_confusion(aldex3_5_var, "ALDEx3", FDR=0.05)
  # limma.a1 <- get_confusion(limma_a1, "limma")
  # limma.a2 <- get_confusion(limma_a2, "limma")
  # limma.a3 <- get_confusion(limma_a3, "limma")
  # limma.90 <- get_confusion(lim.90, "limma")
  # limma.80 <- get_confusion(lim.80, "limma")
  # limma.70 <- get_confusion(lim.70, "limma")
  
  # coeff <- as.numeric(rownames(ald0.conf$TPFPR))
  # par(mfrow=c(1,2))
  
  # take the get_confusion output and plot all TPR/FDR lines
  
  # ald0.conf is the output of get_confusion
  coef <- as.numeric(names(ald0.conf$diff_zero))
  
  TPR_raw <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_raw <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  TPR_0 <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_0 <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  
  TPR_raw3 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_raw3 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  TPR_03 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_03 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  
  for(i in 1:length(coef)){
    TPR_raw[,i] <- ald0.conf$raw_coeff[[i]]$TPR
    FDR_raw[,i] <- ald0.conf$raw_coeff[[i]]$FDR
    TPR_0[,i] <- ald0.conf$raw_zero[[i]]$TPR
    FDR_0[,i] <- ald0.conf$raw_zero[[i]]$FDR
    
    TPR_raw3[,i] <- ald3_0.conf$raw_coeff[[i]]$TPR
    FDR_raw3[,i] <- ald3_0.conf$raw_coeff[[i]]$FDR
    TPR_03[,i] <- ald3_0.conf$raw_zero[[i]]$TPR
    FDR_03[,i] <- ald3_0.conf$raw_zero[[i]]$FDR
  }
  
  par(mfrow=c(1,1))
  
  plot(coef, TPR_raw[1,], col=rgb(0,0,0,0.2), ylim=c(0,1), type='l')
  points(coef, FDR_raw[1,], col=rgb(0,0,1,0.2), type='l')
  points(coef, TPR_0[1,], col=rgb(1,0,0,0.2), type='l')
  points(coef, FDR_0[1,], col=rgb(1,0,1,0.2), type='l')
  points(coef, TPR_raw3[1,], col=rgb(0,0.2,0.4,0.2))
  points(coef, FDR_raw3[1,], col=rgb(0,0.6,1,0.2), type='l')
  points(coef, TPR_03[1,], col=rgb(1,0.2,0.4,0.2), type='l')
  points(coef, FDR_03[1,], col=rgb(0.8,0.2,0.6,0.2), type='l')
  
  for(i in 2:nrow(TPR_raw)){
    points(coef, TPR_raw[i,], col=rgb(0,0,0,0.2), type='l')
    points(coef, FDR_raw[i,], col=rgb(0,0,1,0.2), type='l')
    points(coef, TPR_0[i,], col=rgb(1,0,0,0.2), type='l')
    points(coef, FDR_0[i,], col=rgb(1,0,1,0.2), type='l')
    points(coef, TPR_raw3[i,], col=rgb(0,0.2,0.4,0.2), type='l')
    points(coef, FDR_raw3[i,], col=rgb(0,0.6,1,0.2), type='l')
    points(coef, TPR_03[i,], col=rgb(1,0.2,0.4,0.2), type='l')
    points(coef, FDR_03[i,], col=rgb(0.8,0.2,0.6,0.2), type='l')
  }
  

#   #plot(coeff, ald0.conf$TPFPR[1,1], col="black", type='b', lty=1, 
#        xlab="modeled LFC", ylab="TPR/FDR", xlim=c(0,1), ylim=c(0,1),
#        main="real LFC > coeff")
#   # points(coeff,ald1.conf$TPFPR[,1], col="black", type='b', pch='1', lty=1, cex=0.6)
#   # points(coeff,ald2.conf$TPFPR[,1], col="black", type='b', pch='2', lty=1, cex=0.6)
#   # points(coeff,ald5.conf$TPFPR[,1], col="black", type='b', pch='5', lty=1, cex=0.6)
#   # points(coeff,ald0.conf$TPFPR[,2], col="black", type='b', lty=2)
#   # points(coeff,ald1.conf$TPFPR[,2], col="black", type='b', lty=2, pch='1', cex=0.6)
#   # points(coeff,ald2.conf$TPFPR[,2], col="black", type='b', lty=2, pch='2', cex=0.6)
#   # points(coeff,ald5.conf$TPFPR[,2], col="black", type='b', lty=2, pch='5', cex=0.6)
#   # 
#   # 
#   # points(coeff, ald0.conf5$TPFPR[,1], col="blue", type='b', lty=1)
#   # points(coeff,ald1.conf5$TPFPR[,1], col="blue", type='b', pch='1', lty=1, cex=0.6)
#   # points(coeff,ald2.conf5$TPFPR[,1], col="blue", type='b', pch='2', lty=1, cex=0.6)
#   # points(coeff,ald5.conf5$TPFPR[,1], col="blue", type='b', pch='5', lty=1, cex=0.6)
#   # points(coeff,ald0.conf5$TPFPR[,2], col="blue", type='b', lty=2)
#   # points(coeff,ald1.conf5$TPFPR[,2], col="blue", type='b', lty=2, pch='1', cex=0.6)
#   # points(coeff,ald2.conf5$TPFPR[,2], col="blue", type='b', lty=2, pch='2', cex=0.6)
#   # points(coeff,ald5.conf5$TPFPR[,2], col="blue", type='b', lty=2, pch='5', cex=0.6)
#   # 
#   # # points(coeff,ald0.conf$TPFPR[,3], col="grey", type='b', lty=2)
#   # # points(coeff,ald1.conf$TPFPR[,3], col="grey", type='b', lty=2, pch='1', cex=0.6)
#   # # points(coeff,ald2.conf$TPFPR[,3], col="grey", type='b', lty=2, pch='2', cex=0.6)
#   # # points(coeff,ald5.conf$TPFPR[,3], col="grey", type='b', lty=2, pch='5', cex=0.6)
#   # # points(coeff,ald0.conf$TPFPR[,4], col="grey", type='b', lty=2)
#   # # points(coeff,ald1.conf$TPFPR[,4], col="grey", type='b', lty=2, pch='1', cex=0.6)
#   # points(coeff,ald2.conf$TPFPR[,4], col="grey", type='b', lty=2, pch='2', cex=0.6)
#   # # points(coeff,ald5.conf$TPFPR[,4], col="grey", type='b', lty=2, pch='5', cex=0.6)
#   # 
#   # points(coeff, ald3_0.conf$TPFPR[,1], col="darkgreen", type='b', lty=1)
#   # points(coeff,ald3_1.conf$TPFPR[,1], col="darkgreen", type='b', pch='1', lty=1, cex=0.6)
#   # points(coeff,ald3_2.conf$TPFPR[,1], col="darkgreen", type='b', pch='2', lty=1, cex=0.6)
#   # points(coeff,ald3_5.conf$TPFPR[,1], col="darkgreen", type='b', pch='5', lty=1, cex=0.6)
#   # points(coeff,ald3_0.conf$TPFPR[,2], col="darkgreen", type='b', lty=2)
#   # points(coeff,ald3_1.conf$TPFPR[,2], col="darkgreen", type='b', lty=2, pch='1', cex=0.6)
#   # points(coeff,ald3_2.conf$TPFPR[,2], col="darkgreen", type='b', lty=2, pch='2', cex=0.6)
#   # points(coeff,ald3_5.conf$TPFPR[,2], col="darkgreen", type='b', lty=2, pch='5', cex=0.6)
#   # 
#   # points(coeff, ald3_0.conf5$TPFPR[,1], col="red", type='b', lty=1)
#   # points(coeff,ald3_1.conf5$TPFPR[,1], col="red", type='b', pch='1', lty=1, cex=0.6)
#   # points(coeff,ald3_2.conf5$TPFPR[,1], col="red", type='b', pch='2', lty=1, cex=0.6)
#   # points(coeff,ald3_5.conf5$TPFPR[,1], col="red", type='b', pch='5', lty=1, cex=0.6)
#   # points(coeff,ald3_0.conf5$TPFPR[,2], col="red", type='b', lty=2)
#   # points(coeff,ald3_1.conf5$TPFPR[,2], col="red", type='b', lty=2, pch='1', cex=0.6)
#   # points(coeff,ald3_2.conf5$TPFPR[,2], col="red", type='b', lty=2, pch='2', cex=0.6)
#   # points(coeff,ald3_5.conf5$TPFPR[,2], col="red", type='b', lty=2, pch='5', cex=0.6)
#   # 
#   # 
#   # # points(coeff,ald3_0.conf$TPFPR[,3], col="lightgreen", type='b', lty=1)
#   # # points(coeff,ald3_1.conf$TPFPR[,3], col="lightgreen", type='b', lty=2, pch='1', cex=0.6)
#   # # points(coeff,ald3_2.conf$TPFPR[,3], col="lightgreen", type='b', lty=2, pch='2', cex=0.6)
#   # # points(coeff,ald3_5.conf$TPFPR[,3], col="lightgreen", type='b', lty=2, pch='5', cex=0.6)
#   # # points(coeff,ald3_0.conf$TPFPR[,4], col="lightgreen", type='b', lty=2)
#   # # points(coeff,ald3_1.conf$TPFPR[,4], col="lightgreen", type='b', lty=2, pch='1', cex=0.6)
#   # # points(coeff,ald3_2.conf$TPFPR[,4], col="lightgreen", type='b', lty=2, pch='2', cex=0.6)
#   # # points(coeff,ald3_5.conf$TPFPR[,4], col="lightgreen", type='b', lty=2, pch='5', cex=0.6)
#   # # 
#   # 
#   # # points(coeff, des.conf$TPFPR[,1], col="blue", type='b', lty=1)
#   # # points(coeff, des.conf$TPFPR[,2], col="blue", type='b', lty=2)
#   # # points(coeff, des.conf$TPFPR[,3], col="steelblue1", type='b', lty=1)
#   # # points(coeff, des.conf$TPFPR[,4], col="steelblue1", type='b', lty=2)
#   # # points(coeff, edge.conf$TPFPR[,1], col="violetred", type='b', lty=1)
#   # # points(coeff, edge.conf$TPFPR[,2], col="violetred", type='b', lty=2)
#   # # points(coeff, edge.conf$TPFPR[,3], col="pink", type='b', lty=1)
#   # # points(coeff, edge.conf$TPFPR[,4], col="pink", type='b', lty=2)
#   #  points(coeff, limma.conf$TPFPR[,1], col="brown", type='b', lty=1)
#   #  points(coeff, limma.conf$TPFPR[,2], col="brown", type='b', lty=2)
#   #  points(coeff, limma.conf5$TPFPR[,1], col="orange", type='b', lty=1)
#   #  points(coeff, limma.conf5$TPFPR[,2], col="orange", type='b', lty=2)
#   # # points(coeff, limma.conf$TPFPR[,3], col="orange", type='b', lty=1)
#   # # points(coeff, limma.conf$TPFPR[,4], col="orange", type='b', lty=2)
#   # 
#   # # points(coeff, limma.a1$TPFPR[,1], col="darkgreen", type='b', lty=1)
#   # # points(coeff, limma.a1$TPFPR[,2], col="darkgreen", type='b', lty=2)
#   # # points(coeff, limma.a1$TPFPR[,3], col="lightgreen", type='b', lty=1)
#   # # points(coeff, limma.a1$TPFPR[,4], col="lightgreen", type='b', lty=2)
#   # # 
#   # # points(coeff, limma.a2$TPFPR[,1], col="yellow2", type='b', lty=1)
#   # # points(coeff, limma.a2$TPFPR[,2], col="yellow2", type='b', lty=2)
#   # points(coeff, limma.a2$TPFPR[,3], col="yellow", type='b', lty=1)
#   # points(coeff, limma.a2$TPFPR[,4], col="yellow", type='b', lty=2)
#   # 
#   # points(coeff, limma.a3$TPFPR[,1], col="purple", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,2], col="purple", type='b', lty=2)
#   # points(coeff, limma.a3$TPFPR[,3], col="violet", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,4], col="violet", type='b', lty=2)
#   # 
#   # points(coeff, lim.90$TPFPR[,1], col= "coral3", type='b', lty=1)
#   # points(coeff, lim.90$TPFPR[,2], col="coral3", type='b', lty=2)
#   # points(coeff, lim.90$TPFPR[,3], col="coral", type='b', lty=1)
#   # points(coeff, lim.90$TPFPR[,4], col="coral", type='b', lty=2)
#   # 
#   # points(coeff, lim.80$TPFPR[,1], col="hotpink2", type='b', lty=1)
#   # points(coeff, lim.80$TPFPR[,2], col="hotpink2", type='b', lty=2)
#   # points(coeff, lim.80$TPFPR[,3], col="hotpink", type='b', lty=1)
#   # points(coeff, lim.80$TPFPR[,4], col= "hotpink", type='b', lty=2)
#   # 
#   # points(coeff, limma.a3$TPFPR[,1], col="cyan3", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,2], col="cyan3", type='b', lty=2)
#   # points(coeff, limma.a3$TPFPR[,3], col="cyan", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,4], col="cyan", type='b', lty=2)
# 
#   abline(h=0.05, lty=3, lwd=2, col="grey") 
#   
#   plot(coeff, ald0.conf$TPFPR[,5], col="black", type='b', lty=1, 
#        xlab="modeled LFC", ylab="TPR/FDR", xlim=c(0,1), ylim=c(0,1),
#        main="real LFC > 0")
# #   points(coeff,ald1.conf$TPFPR[,5], col="black", type='b', pch='1', lty=1, cex=0.6)
# #   points(coeff,ald2.conf$TPFPR[,5], col="black", type='b', pch='2', lty=1, cex=0.6)
# #   points(coeff,ald5.conf$TPFPR[,5], col="black", type='b', pch='5', lty=1, cex=0.6)
# #   
# #   points(coeff,ald0.conf$TPFPR[,6], col="black", type='b', lty=2)
# #   points(coeff,ald1.conf$TPFPR[,6], col="black", type='b', lty=2, pch='1', cex=0.6)
# #   points(coeff,ald2.conf$TPFPR[,6], col="black", type='b', lty=2, pch='2', cex=0.6)
# #   points(coeff,ald5.conf$TPFPR[,6], col="black", type='b', lty=2, pch='5', cex=0.6)
# #   
# #   points(coeff, ald0.conf5$TPFPR[,5], col="blue", type='b', lty=1)
# #   points(coeff,ald1.conf5$TPFPR[,5], col="blue", type='b', pch='1', lty=1, cex=0.6)
# #   points(coeff,ald2.conf5$TPFPR[,5], col="blue", type='b', pch='2', lty=1, cex=0.6)
# #   points(coeff,ald5.conf5$TPFPR[,5], col="blue", type='b', pch='5', lty=1, cex=0.6)
# #   points(coeff,ald0.conf5$TPFPR[,6], col="blue", type='b', lty=2)
# #   points(coeff,ald1.conf5$TPFPR[,6], col="blue", type='b', lty=2, pch='1', cex=0.6)
# #   points(coeff,ald2.conf5$TPFPR[,6], col="blue", type='b', lty=2, pch='2', cex=0.6)
# #   points(coeff,ald5.conf5$TPFPR[,6], col="blue", type='b', lty=2, pch='5', cex=0.6)
# #   
# #   # points(coeff,ald0.conf$TPFPR[,7], col="grey", type='b', lty=2)
# #   # points(coeff,ald1.conf$TPFPR[,7], col="grey", type='b', lty=2, pch='1', cex=0.6)
# #   # points(coeff,ald2.conf$TPFPR[,7], col="grey", type='b', lty=2, pch='2', cex=0.6)
# #   # points(coeff,ald5.conf$TPFPR[,7], col="grey", type='b', lty=2, pch='5', cex=0.6)
# #   # 
# #   # points(coeff,ald0.conf$TPFPR[,8], col="grey", type='b', lty=2)
# #   # points(coeff,ald1.conf$TPFPR[,8], col="grey", type='b', lty=2, pch='1', cex=0.6)
# #   # points(coeff,ald2.conf$TPFPR[,8], col="grey", type='b', lty=2, pch='2', cex=0.6)
# #   # points(coeff,ald5.conf$TPFPR[,8], col="grey", type='b', lty=2, pch='5', cex=0.6)
# #   # abline(h=0.05, lty=3, lwd=2, col="grey")
# #   
# #   points(coeff, ald3_0.conf$TPFPR[,5], col="darkgreen", type='b', lty=2)
# #   points(coeff,ald3_1.conf$TPFPR[,5], col="darkgreen", type='b', pch='1', lty=1, cex=0.6)
# #   points(coeff,ald3_2.conf$TPFPR[,5], col="darkgreen", type='b', pch='2', lty=1, cex=0.6)
# #   points(coeff,ald3_5.conf$TPFPR[,5], col="darkgreen", type='b', pch='5', lty=1, cex=0.6)
# #   points(coeff,ald3_0.conf$TPFPR[,6], col="darkgreen", type='b', lty=2)
# #   points(coeff,ald3_1.conf$TPFPR[,6], col="darkgreen", type='b', lty=2, pch='1', cex=0.6)
# #   points(coeff,ald3_2.conf$TPFPR[,6], col="darkgreen", type='b', lty=2, pch='2', cex=0.6)
# #   points(coeff,ald3_5.conf$TPFPR[,6], col="darkgreen", type='b', lty=2, pch='5', cex=0.6)
# #   
# #   points(coeff, ald3_0.conf5$TPFPR[,5], col="red", type='b', lty=1)
# #   points(coeff,ald3_1.conf5$TPFPR[,5], col="red", type='b', pch='1', lty=1, cex=0.6)
# #   points(coeff,ald3_2.conf5$TPFPR[,5], col="red", type='b', pch='2', lty=1, cex=0.6)
# #   points(coeff,ald3_5.conf5$TPFPR[,5], col="red", type='b', pch='5', lty=1, cex=0.6)
# #   points(coeff,ald3_0.conf5$TPFPR[,6], col="red", type='b', lty=2)
# #   points(coeff,ald3_1.conf5$TPFPR[,6], col="red", type='b', lty=2, pch='1', cex=0.6)
# #   points(coeff,ald3_2.conf5$TPFPR[,6], col="red", type='b', lty=2, pch='2', cex=0.6)
# #   points(coeff,ald3_5.conf5$TPFPR[,6], col="red", type='b', lty=2, pch='5', cex=0.6)
# #   
# #   # points(coeff,ald3_0.conf$TPFPR[,7], col="lightgreen", type='b', lty=2)
# #   # points(coeff,ald3_1.conf$TPFPR[,7], col="lightgreen", type='b', lty=2, pch='1', cex=0.6)
# #   # points(coeff,ald3_2.conf$TPFPR[,7], col="lightgreen", type='b', lty=2, pch='2', cex=0.6)
# #   # points(coeff,ald3_5.conf$TPFPR[,7], col="lightgreen", type='b', lty=2, pch='5', cex=0.6)
# #   # points(coeff,ald3_0.conf$TPFPR[,8], col="lightgreen", type='b', lty=2)
# #   # points(coeff,ald3_1.conf$TPFPR[,8], col="lightgreen", type='b', lty=2, pch='1', cex=0.6)
# #   # points(coeff,ald3_2.conf$TPFPR[,8], col="lightgreen", type='b', lty=2, pch='2', cex=0.6)
# #   # points(coeff,ald3_5.conf$TPFPR[,8], col="lightgreen", type='b', lty=2, pch='5', cex=0.6)
# # # 
# # #   points(coeff, des.conf$TPFPR[,5], col="blue", type='b', lty=1)
# # #   points(coeff, des.conf$TPFPR[,6], col="blue", type='b', lty=2)
# # #   points(coeff, des.conf$TPFPR[,7], col="steelblue1", type='b', lty=1)
# # #   points(coeff, des.conf$TPFPR[,8], col="steelblue1", type='b', lty=2)
# # #   points(coeff, edge.conf$TPFPR[,5], col="violetred", type='b', lty=1)
# # #   points(coeff, edge.conf$TPFPR[,6], col="violetred", type='b', lty=2)
# # #   points(coeff, edge.conf$TPFPR[,7], col="pink", type='b', lty=1)
# # #   points(coeff, edge.conf$TPFPR[,8], col="pink", type='b', lty=2)
# #    points(coeff, limma.conf$TPFPR[,5], col="brown", type='b', lty=1)
# #    points(coeff, limma.conf$TPFPR[,6], col="brown", type='b', lty=2)
# #    points(coeff, limma.conf5$TPFPR[,5], col="orange", type='b', lty=1)
# #    points(coeff, limma.conf5$TPFPR[,6], col="orange", type='b', lty=2)
# # #   points(coeff, limma.conf$TPFPR[,7], col="orange", type='b', lty=1)
# # #   points(coeff, limma.conf$TPFPR[,8], col="orange", type='b', lty=2)
# 
#   # points(coeff, limma.a1$TPFPR[,5], col="darkgreen", type='b', lty=1)
#   # points(coeff, limma.a1$TPFPR[,6], col="darkgreen", type='b', lty=2)
#   # points(coeff, limma.a1$TPFPR[,7], col="lightgreen", type='b', lty=1)
#   # points(coeff, limma.a1$TPFPR[,8], col="lightgreen", type='b', lty=2)
#   # 
#   # points(coeff, limma.a2$TPFPR[,1], col="yellow2", type='b', lty=1)
#   # points(coeff, limma.a2$TPFPR[,2], col="yellow2", type='b', lty=2)
#   # points(coeff, limma.a2$TPFPR[,3], col="yellow", type='b', lty=1)
#   # points(coeff, limma.a2$TPFPR[,4], col="yellow", type='b', lty=2)
#   # 
#   # points(coeff, limma.a3$TPFPR[,5], col="purple", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,6], col="purple", type='b', lty=2)
#   # points(coeff, limma.a3$TPFPR[,7], col="violet", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,8], col="violet", type='b', lty=2)
#   # 
#   # points(coeff, lim.90$TPFPR[,5], col= "coral3", type='b', lty=1)
#   # points(coeff, lim.90$TPFPR[,6], col="coral3", type='b', lty=2)
#   # points(coeff, lim.90$TPFPR[,7], col="coral", type='b', lty=1)
#   # points(coeff, lim.90$TPFPR[,8], col="coral", type='b', lty=2)
#   # 
#   # points(coeff, lim.80$TPFPR[,5], col="hotpink2", type='b', lty=1)
#   # points(coeff, lim.80$TPFPR[,6], col="hotpink2", type='b', lty=2)
#   # points(coeff, lim.80$TPFPR[,7], col="hotpink", type='b', lty=1)
#   # points(coeff, lim.80$TPFPR[,8], col= "hotpink", type='b', lty=2)
#   # 
#   # points(coeff, limma.a3$TPFPR[,5], col="cyan3", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,6], col="cyan3", type='b', lty=2)
#   # points(coeff, limma.a3$TPFPR[,7], col="cyan", type='b', lty=1)
#   # points(coeff, limma.a3$TPFPR[,8], col="cyan", type='b', lty=2)
#   # 
#   
#   abline(h=0.05, lty=3, lwd=2, col="grey")
#   
# }
  
