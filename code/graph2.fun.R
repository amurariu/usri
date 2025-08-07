#graphing function for panel containing 0.5 log fold change cutoff

source('code/get_confusion.R')
#DESeq_var, edgeR_var, limma_var,

graph2.fun <- function(DESeq_var, edgeR_var, limma_var, aldex_0_var, aldex_1_var, aldex_2_var, aldex_5_var, aldex3_0_var, aldex3_1_var, aldex3_2_var, aldex3_5_var){
  
#perform get_confusion function first
des.conf <- get_confusion(DESeq_var, "DESeq", FDR=0.05)
edge.conf <- get_confusion(edgeR_var, "edgeR", FDR=0.05)
ald0.conf <- get_confusion(aldex_0_var, "ALDEx2", FDR=0.05)
ald1.conf <- get_confusion(aldex_1_var, "ALDEx2", FDR=0.05)
ald2.conf <- get_confusion(aldex_2_var, "ALDEx2", FDR=0.05)
ald5.conf <- get_confusion(aldex_5_var, "ALDEx2", FDR=0.05)
ald3_0.conf <- get_confusion(aldex3_0_var, "ALDEx3", FDR=0.05)
ald3_1.conf <- get_confusion(aldex3_1_var, "ALDEx3", FDR=0.05)
ald3_2.conf <- get_confusion(aldex3_2_var, "ALDEx3", FDR=0.05)
ald3_5.conf <- get_confusion(aldex3_5_var, "ALDEx3", FDR=0.05)
limma.conf <- get_confusion(limma_var, "limma", FDR=0.05)
  
coeff <- as.numeric(rownames(ald0.conf$TPFPR))
par(mfrow=c(2,2))
  
#left-side panel

  plot(coeff, ald0.conf$TPFPR[,1], col="black", type='b', lty=1,
       xlab="modeled LFC", ylab="TPR", xlim=c(0,1), ylim=c(0,1),
       main="TPR with real LFC > coeff")
  points(coeff,ald1.conf$TPFPR[,1], col="black", type='b', pch='1', lty=1, cex=0.6) #plotting cTPR0 from TPFPR
  points(coeff,ald2.conf$TPFPR[,1], col="black", type='b', pch='2', lty=1, cex=0.6)
  points(coeff,ald5.conf$TPFPR[,1], col="black", type='b', pch='5', lty=1, cex=0.6)
  # points(coeff,ald0.conf$TPFPR[,2], col="black", type='b', lty=1) #plotting cFDR0 from TPFPR
  # points(coeff,ald1.conf$TPFPR[,2], col="black", type='b', lty=1, pch='1', cex=0.6) 
  # points(coeff,ald2.conf$TPFPR[,2], col="black", type='b', lty=1, pch='2', cex=0.6)
  # points(coeff,ald5.conf$TPFPR[,2], col="black", type='b', lty=1, pch='5', cex=0.6)
  
  points(coeff,ald0.conf$TPFPR[,3], col="grey", type='b', lty=3) #plotting cTPR5 from TPFPR
  points(coeff,ald1.conf$TPFPR[,3], col="grey", type='b', lty=3, pch='1', cex=0.6)
  points(coeff,ald2.conf$TPFPR[,3], col="grey", type='b', lty=3, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,3], col="grey", type='b', lty=3, pch='5', cex=0.6)
  # points(coeff,ald0.conf$TPFPR[,4], col="grey", type='b', lty=3) #plotting cFDR5 from TPFPR
  # points(coeff,ald1.conf$TPFPR[,4], col="grey", type='b', lty=3, pch='1', cex=0.6)
  # points(coeff,ald2.conf$TPFPR[,4], col="grey", type='b', lty=3, pch='2', cex=0.6)
  # points(coeff,ald5.conf$TPFPR[,4], col="grey", type='b', lty=3, pch='5', cex=0.6)
  abline(h=0.05, lty=3, lwd=2, col="grey")
  
  points(coeff, ald3_0.conf$TPFPR[,1], col="blue", type='b', lty=1) #plotting cTPR5 from TPFPR
  points(coeff,ald3_1.conf$TPFPR[,1], col="blue", type='b', pch='1', lty=1, cex=0.6)
  points(coeff,ald3_2.conf$TPFPR[,1], col="blue", type='b', pch='2', lty=1, cex=0.6)
  points(coeff,ald3_5.conf$TPFPR[,1], col="blue", type='b', pch='5', lty=1, cex=0.6)
  # points(coeff,ald3_0.conf$TPFPR[,2], col="blue", type='b', lty=1) #plotting cFDR5 from TPFPR
  # points(coeff,ald3_1.conf$TPFPR[,2], col="blue", type='b', lty=1, pch='1', cex=0.6)
  # points(coeff,ald3_2.conf$TPFPR[,2], col="blue", type='b', lty=1, pch='2', cex=0.6)
  # points(coeff,ald3_5.conf$TPFPR[,2], col="blue", type='b', lty=1, pch='5', cex=0.6)
  points(coeff,ald3_0.conf$TPFPR[,3], col="lightblue", type='b', lty=3)
  points(coeff,ald3_1.conf$TPFPR[,3], col="lightblue", type='b', lty=3, pch='1', cex=0.6)
  points(coeff,ald3_2.conf$TPFPR[,3], col="lightblue", type='b', lty=3, pch='2', cex=0.6)
  points(coeff,ald3_5.conf$TPFPR[,3], col="lightblue", type='b', lty=3, pch='5', cex=0.6)
  # points(coeff,ald3_0.conf$TPFPR[,4], col="steelblue1", type='b', lty=3)
  # points(coeff,ald3_1.conf$TPFPR[,4], col="steelblue1", type='b', lty=3, pch='1', cex=0.6)
  # points(coeff,ald3_2.conf$TPFPR[,4], col="steelblue1", type='b', lty=3, pch='2', cex=0.6)
  # points(coeff,ald3_5.conf$TPFPR[,4], col="steelblue1", type='b', lty=3, pch='5', cex=0.6)
  # 
  points(coeff, des.conf$TPFPR[,1], col="gold3", type='b', lty=1)
  # points(coeff, des.conf$TPFPR[,2], col="gold3", type='b', lty=1)
  points(coeff, des.conf$TPFPR[,3], col="gold1", type='b', lty=3)
  # points(coeff, des.conf$TPFPR[,4], col="gold1", type='b', lty=3)
  points(coeff, edge.conf$TPFPR[,1], col="hotpink", type='b', lty=1)
  # points(coeff, edge.conf$TPFPR[,2], col="violetred", type='b', lty=1)
  points(coeff, edge.conf$TPFPR[,3], col="pink", type='b', lty=3)
  # points(coeff, edge.conf$TPFPR[,4], col="pink", type='b', lty=3)
  points(coeff, limma.conf$TPFPR[,1], col="brown", type='b', lty=1)
  # points(coeff, limma.conf$TPFPR[,2], col="brown", type='b', lty=1)
  points(coeff, limma.conf$TPFPR[,3], col="orange", type='b', lty=3)
  # points(coeff, limma.conf$TPFPR[,4], col="orange", type='b', lty=3)
  abline(h=0.05, lty=3, lwd=2, col="grey")

  
#right side panel
plot(coeff, ald0.conf$TPFPR[,5], col="black", type='b', lty=1,
     xlab="modeled LFC", ylab="TPR", xlim=c(0,1), ylim=c(0,1),
     main="TPR with real LFC > 0")
points(coeff,ald1.conf$TPFPR[,5], col="black", type='b', pch='1', lty=1, cex=0.6) #plotting cTPR0 from TPFPR
points(coeff,ald2.conf$TPFPR[,5], col="black", type='b', pch='2', lty=1, cex=0.6)
points(coeff,ald5.conf$TPFPR[,5], col="black", type='b', pch='5', lty=1, cex=0.6)
# points(coeff,ald0.conf$TPFPR[,6], col="black", type='b', lty=1) #plotting cFDR0 from TPFPR
# points(coeff,ald1.conf$TPFPR[,6], col="black", type='b', lty=1, pch='1', cex=0.6) 
# points(coeff,ald2.conf$TPFPR[,6], col="black", type='b', lty=1, pch='2', cex=0.6)
# points(coeff,ald5.conf$TPFPR[,6], col="black", type='b', lty=1, pch='5', cex=0.6)

points(coeff,ald0.conf$TPFPR[,7], col="grey", type='b', lty=3) #plotting cTPR5 from TPFPR
points(coeff,ald1.conf$TPFPR[,7], col="grey", type='b', lty=3, pch='1', cex=0.6)
points(coeff,ald2.conf$TPFPR[,7], col="grey", type='b', lty=3, pch='2', cex=0.6)
points(coeff,ald5.conf$TPFPR[,7], col="grey", type='b', lty=3, pch='5', cex=0.6)
# points(coeff,ald0.conf$TPFPR[,8], col="grey", type='b', lty=3) #plotting cFDR5 from TPFPR
# points(coeff,ald1.conf$TPFPR[,8], col="grey", type='b', lty=3, pch='1', cex=0.6)
# points(coeff,ald2.conf$TPFPR[,8], col="grey", type='b', lty=3, pch='2', cex=0.6)
# points(coeff,ald5.conf$TPFPR[,8], col="grey", type='b', lty=3, pch='5', cex=0.6)
abline(h=0.05, lty=3, lwd=2, col="grey")

points(coeff, ald3_0.conf$TPFPR[,5], col="blue", type='b', lty=1) #plotting cTPR5 from TPFPR
points(coeff,ald3_1.conf$TPFPR[,5], col="blue", type='b', pch='1', lty=1, cex=0.6)
points(coeff,ald3_2.conf$TPFPR[,5], col="blue", type='b', pch='2', lty=1, cex=0.6)
points(coeff,ald3_5.conf$TPFPR[,5], col="blue", type='b', pch='5', lty=1, cex=0.6)
# points(coeff,ald3_0.conf$TPFPR[,6], col="blue", type='b', lty=1) #plotting cFDR5 from TPFPR
# points(coeff,ald3_1.conf$TPFPR[,6], col="blue", type='b', lty=1, pch='1', cex=0.6)
# points(coeff,ald3_2.conf$TPFPR[,6], col="blue", type='b', lty=1, pch='2', cex=0.6)
# points(coeff,ald3_5.conf$TPFPR[,6], col="blue", type='b', lty=1, pch='5', cex=0.6)
points(coeff,ald3_0.conf$TPFPR[,7], col="lightblue", type='b', lty=3)
points(coeff,ald3_1.conf$TPFPR[,7], col="lightblue", type='b', lty=3, pch='1', cex=0.6)
points(coeff,ald3_2.conf$TPFPR[,7], col="lightblue", type='b', lty=3, pch='2', cex=0.6)
points(coeff,ald3_5.conf$TPFPR[,7], col="lightblue", type='b', lty=3, pch='5', cex=0.6)
# points(coeff,ald3_0.conf$TPFPR[,8], col="steelblue1", type='b', lty=3)
# points(coeff,ald3_1.conf$TPFPR[,8], col="steelblue1", type='b', lty=3, pch='1', cex=0.6)
# points(coeff,ald3_2.conf$TPFPR[,8], col="steelblue1", type='b', lty=3, pch='2', cex=0.6)
# points(coeff,ald3_5.conf$TPFPR[,8], col="steelblue1", type='b', lty=3, pch='5', cex=0.6)

points(coeff, des.conf$TPFPR[,5], col="gold3", type='b', lty=1)
# points(coeff, des.conf$TPFPR[,6], col="gold3", type='b', lty=1)
points(coeff, des.conf$TPFPR[,7], col="gold1", type='b', lty=3)
# points(coeff, des.conf$TPFPR[,8], col="gold1", type='b', lty=3)
points(coeff, edge.conf$TPFPR[,5], col="hotpink", type='b', lty=1)
# points(coeff, edge.conf$TPFPR[,6], col="violetred", type='b', lty=1)
points(coeff, edge.conf$TPFPR[,7], col="pink", type='b', lty=3)
# points(coeff, edge.conf$TPFPR[,8], col="pink", type='b', lty=3)
points(coeff, limma.conf$TPFPR[,5], col="brown", type='b', lty=1)
# points(coeff, limma.conf$TPFPR[,6], col="brown", type='b', lty=1)
points(coeff, limma.conf$TPFPR[,7], col="orange", type='b', lty=3)
# points(coeff, limma.conf$TPFPR[,8], col="orange", type='b', lty=3)
abline(h=0.05, lty=3, lwd=2, col="grey")


#bottom left panel
plot(coeff, ald0.conf$TPFPR[,2], col="black", type='b', lty=1,
     xlab="modeled LFC", ylab="FDR", xlim=c(0,1), ylim=c(0,1),
     main="FDR with real LFC > coeff")
points(coeff,ald1.conf$TPFPR[,2], col="black", type='b', lty=1, pch='1', cex=0.6) 
points(coeff,ald2.conf$TPFPR[,2], col="black", type='b', lty=1, pch='2', cex=0.6)
points(coeff,ald5.conf$TPFPR[,2], col="black", type='b', lty=1, pch='5', cex=0.6)

points(coeff,ald0.conf$TPFPR[,4], col="grey", type='b', lty=3) #plotting cFDR5 from TPFPR
points(coeff,ald1.conf$TPFPR[,4], col="grey", type='b', lty=3, pch='1', cex=0.6)
points(coeff,ald2.conf$TPFPR[,4], col="grey", type='b', lty=3, pch='2', cex=0.6)
points(coeff,ald5.conf$TPFPR[,4], col="grey", type='b', lty=3, pch='5', cex=0.6)

points(coeff,ald3_0.conf$TPFPR[,2], col="blue", type='b', lty=1) #plotting cFDR5 from TPFPR
points(coeff,ald3_1.conf$TPFPR[,2], col="blue", type='b', lty=1, pch='1', cex=0.6)
points(coeff,ald3_2.conf$TPFPR[,2], col="blue", type='b', lty=1, pch='2', cex=0.6)
points(coeff,ald3_5.conf$TPFPR[,2], col="blue", type='b', lty=1, pch='5', cex=0.6)

points(coeff,ald3_0.conf$TPFPR[,4], col="lightblue", type='b', lty=3)
points(coeff,ald3_1.conf$TPFPR[,4], col="lightblue", type='b', lty=3, pch='1', cex=0.6)
points(coeff,ald3_2.conf$TPFPR[,4], col="lightblue", type='b', lty=3, pch='2', cex=0.6)
points(coeff,ald3_5.conf$TPFPR[,4], col="lightblue", type='b', lty=3, pch='5', cex=0.6)

points(coeff, des.conf$TPFPR[,2], col="gold3", type='b', lty=1)
points(coeff, des.conf$TPFPR[,4], col="gold1", type='b', lty=3)
points(coeff, edge.conf$TPFPR[,2], col="hotpink", type='b', lty=1)
points(coeff, edge.conf$TPFPR[,4], col="pink", type='b', lty=3)
points(coeff, limma.conf$TPFPR[,2], col="brown", type='b', lty=1)
points(coeff, limma.conf$TPFPR[,4], col="orange", type='b', lty=3)
abline(h=0.05, lty=3, lwd=2, col="grey")

#bottom right panel

#right side panel
plot(coeff, ald0.conf$TPFPR[,6], col="black", type='b', lty=1,
     xlab="modeled LFC", ylab="FDR", xlim=c(0,1), ylim=c(0,1),
     main="FDR with real LFC > 0")
points(coeff,ald1.conf$TPFPR[,6], col="black", type='b', lty=1, pch='1', cex=0.6) 
points(coeff,ald2.conf$TPFPR[,6], col="black", type='b', lty=1, pch='2', cex=0.6)
points(coeff,ald5.conf$TPFPR[,6], col="black", type='b', lty=1, pch='5', cex=0.6)

points(coeff,ald0.conf$TPFPR[,8], col="grey", type='b', lty=3) #plotting cFDR5 from TPFPR
points(coeff,ald1.conf$TPFPR[,8], col="grey", type='b', lty=3, pch='1', cex=0.6)
points(coeff,ald2.conf$TPFPR[,8], col="grey", type='b', lty=3, pch='2', cex=0.6)
points(coeff,ald5.conf$TPFPR[,8], col="grey", type='b', lty=3, pch='5', cex=0.6)

points(coeff,ald3_0.conf$TPFPR[,6], col="blue", type='b', lty=1) #plotting cFDR5 from TPFPR
points(coeff,ald3_1.conf$TPFPR[,6], col="blue", type='b', lty=1, pch='1', cex=0.6)
points(coeff,ald3_2.conf$TPFPR[,6], col="blue", type='b', lty=1, pch='2', cex=0.6)
points(coeff,ald3_5.conf$TPFPR[,6], col="blue", type='b', lty=1, pch='5', cex=0.6)
points(coeff,ald3_0.conf$TPFPR[,8], col="lightblue", type='b', lty=3)
points(coeff,ald3_1.conf$TPFPR[,8], col="lightblue", type='b', lty=3, pch='1', cex=0.6)
points(coeff,ald3_2.conf$TPFPR[,8], col="lightblue", type='b', lty=3, pch='2', cex=0.6)
points(coeff,ald3_5.conf$TPFPR[,8], col="lightblue", type='b', lty=3, pch='5', cex=0.6)

points(coeff, des.conf$TPFPR[,6], col="gold3", type='b', lty=1)
points(coeff, des.conf$TPFPR[,8], col="gold1", type='b', lty=3)
points(coeff, edge.conf$TPFPR[,6], col="hotpink", type='b', lty=1)
points(coeff, edge.conf$TPFPR[,8], col="pink", type='b', lty=3)
points(coeff, limma.conf$TPFPR[,6], col="brown", type='b', lty=1)
points(coeff, limma.conf$TPFPR[,8], col="orange", type='b', lty=3)
abline(h=0.05, lty=3, lwd=2, col="grey")



}



  