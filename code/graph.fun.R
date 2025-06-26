source('code/get_confusion.R')

graph.fun <- function( DESeq_var, edgeR_var, limma_var, aldex_0_var, aldex_1_var, aldex_2_var, aldex_5_var){
  
  #perform get_confusion function first
  des.conf <- get_confusion(DESeq_var, "DESeq")
  edge.conf <- get_confusion(edgeR_var, "edgeR")
  limma.conf <- get_confusion(limma_var, "limma")
  ald0.conf <- get_confusion(aldex_0_var, "ALDEx2")
  ald1.conf <- get_confusion(aldex_1_var, "ALDEx2")
  ald2.conf <- get_confusion(aldex_2_var, "ALDEx2")
  ald5.conf <- get_confusion(aldex_5_var, "ALDEx2")
  
  
  coeff <- as.numeric(rownames(ald0.imm$TPFPR))
  par(mfrow=c(1,2))
  plot(coeff, ald0.conf$TPFPR[,1], col="black", type='b', lty=1, 
       xlab="modeled LFC", ylab="TPR/FDR", xlim=c(0,1), ylim=c(0,1),
       main="FP = min LFC")
  points(coeff,ald1.conf$TPFPR[,1], col="black", type='b', pch='1', lty=1, cex=0.6)
  points(coeff,ald2.conf$TPFPR[,1], col="black", type='b', pch='2', lty=1, cex=0.6)
  points(coeff,ald5.conf$TPFPR[,1], col="black", type='b', pch='5', lty=1, cex=0.6)
  points(coeff,ald0.conf$TPFPR[,2], col="black", type='b', lty=2)
  points(coeff,ald1.conf$TPFPR[,2], col="black", type='b', lty=2, pch='1', cex=0.6))
  points(coeff,ald2.conf$TPFPR[,2], col="black", type='b', lty=2, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,2], col="black", type='b', lty=2, pch='5', cex=0.6)
  
  points(coeff,ald0.conf$TPFPR[,3], col="grey", type='b', lty=2)
  points(coeff,ald1.conf$TPFPR[,3], col="grey", type='b', lty=2, pch='1', cex=0.6)
  points(coeff,ald2.conf$TPFPR[,3], col="grey", type='b', lty=2, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,3], col="grey", type='b', lty=2, pch='5', cex=0.6)
  points(coeff,ald0.conf$TPFPR[,4], col="grey", type='b', lty=2)
  points(coeff,ald1.conf$TPFPR[,4], col="grey", type='b', lty=2, pch='1', cex=0.6)
  points(coeff,ald2.conf$TPFPR[,4], col="grey", type='b', lty=2, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,4], col="grey", type='b', lty=2, pch='5', cex=0.6)
  
  points(coeff, des.conf$TPFPR[,1], col="blue", type='b', lty=1)
  points(coeff, des.conf$TPFPR[,2], col="blue", type='b', lty=2)
  points(coeff, des.conf$TPFPR[,3], col="lightblue", type='b', lty=1)
  points(coeff, des.conf$TPFPR[,4], col="lightblue", type='b', lty=2)
  points(coeff, edge.conf$TPFPR[,1], col="red", type='b', lty=1)
  points(coeff, edge.conf$TPFPR[,2], col="red", type='b', lty=2)
  points(coeff, edge.conf$TPFPR[,3], col="pink", type='b', lty=1)
  points(coeff, edge.conf$TPFPR[,4], col="pink", type='b', lty=2)
  points(coeff, limma.conf$TPFPR[,1], col="brown", type='b', lty=1)
  points(coeff, limma.conf$TPFPR[,2], col="brown", type='b', lty=2)
  points(coeff, limma.conf$TPFPR[,3], col="orange", type='b', lty=1)
  points(coeff, limma.conf$TPFPR[,4], col="orange", type='b', lty=2)
  abline(h=0.05, lty=3, lwd=2, col="grey") ################################################
  
  plot(coeff, ald0.conf$TPFPR[,5], col="black", type='b', lty=1, 
       xlab="modeled LFC", ylab="TPR/FDR", xlim=c(0,1), ylim=c(0,1),
       main="FP = min LFC")
  points(coeff,ald1.conf$TPFPR[,5], col="black", type='b', pch='1', lty=1, cex=0.6)
  points(coeff,ald2.conf$TPFPR[,5], col="black", type='b', pch='2', lty=1, cex=0.6)
  points(coeff,ald5.conf$TPFPR[,5], col="black", type='b', pch='5', lty=1, cex=0.6)
  
  points(coeff,ald0.conf$TPFPR[,6], col="black", type='b', lty=2)
  points(coeff,ald1.conf$TPFPR[,6], col="black", type='b', lty=2, pch='1', cex=0.6)
  points(coeff,ald2.conf$TPFPR[,6], col="black", type='b', lty=2, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,6], col="black", type='b', lty=2, pch='5', cex=0.6)
  
  points(coeff,ald0.conf$TPFPR[,7], col="grey", type='b', lty=2)
  points(coeff,ald1.conf$TPFPR[,7], col="grey", type='b', lty=2, pch='1', cex=0.6)
  points(coeff,ald2.conf$TPFPR[,7], col="grey", type='b', lty=2, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,7], col="grey", type='b', lty=2, pch='5', cex=0.6)
  
  points(coeff,ald0.conf$TPFPR[,8], col="grey", type='b', lty=2)
  points(coeff,ald1.conf$TPFPR[,8], col="grey", type='b', lty=2, pch='1', cex=0.6)
  points(coeff,ald2.conf$TPFPR[,8], col="grey", type='b', lty=2, pch='2', cex=0.6)
  points(coeff,ald5.conf$TPFPR[,8], col="grey", type='b', lty=2, pch='5', cex=0.6)
  abline(h=0.05, lty=3, lwd=2, col="grey")
  
  points(coeff, des.conf$TPFPR[,5], col="blue", type='b', lty=1)
  points(coeff, des.conf$TPFPR[,6], col="blue", type='b', lty=2)
  points(coeff, des.conf$TPFPR[,7], col="lightblue", type='b', lty=1)
  points(coeff, des.conf$TPFPR[,8], col="lightblue", type='b', lty=2)
  points(coeff, edge.conf$TPFPR[,5], col="red", type='b', lty=1)
  points(coeff, edge.conf$TPFPR[,6], col="red", type='b', lty=2)
  points(coeff, edge.conf$TPFPR[,7], col="pink", type='b', lty=1)
  points(coeff, edge.conf$TPFPR[,8], col="pink", type='b', lty=2)
  points(coeff, limma.conf$TPFPR[,5], col="brown", type='b', lty=1)
  points(coeff, limma.conf$TPFPR[,6], col="brown", type='b', lty=2)
  points(coeff, limma.conf$TPFPR[,7], col="orange", type='b', lty=1)
  points(coeff, limma.conf$TPFPR[,8], col="orange", type='b', lty=2)
  abline(h=0.05, lty=3, lwd=2, col="grey")
  
}
  
