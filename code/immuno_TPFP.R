analysis.fun <- function(input=NULL, type=NULL, nloop=1) {

  if (type == "DESeq") {
    padj = 6
  } else if (type == "edgeR") {
    padj = 5
  } else if (type == "aldex0.we") {
    padj = 10 #for welches
  } else if (type == "aldex0.wi") {
    padj = 12  #for wilcoxon
  } else if (type == "aldex0.2.we") {
    padj = 10 
  } else if (type == "aldex0.2.wi") {
    padj = 12 
  } else if (type == "aldex0.5.we") {
    padj = 10 
  } else if (type == "aldex0.5.wi") {
    padj = 12
  }

  
#if (input == NULL) {stop("input modelled data")} 


### wilcoxon (t-test is labled welches below)
# there is duplication for deseq with the we and wilcoxon
#analysis.wi <- matrix(data=NA, nrow=300, ncol=7) #change nrow/ncol?
#analysis.wi <- as.data.frame(analysis.wi)
#colnames(analysis.wi) <- c("coeff","iter", "met","PPV","FDR","SEN","SPE")
#met = c("ald", "ald2", "ald5", "des", "des0.5", "des1")
#analysis.wi[,3] <- rep(met,60)

#### we.eBH
#analysis.we <- matrix(data=NA, nrow=300, ncol=7)
#analysis.we <- as.data.frame(analysis.we)
#colnames(analysis.we) <- c("coeff","iter", "met","PPV","FDR","SEN","SPE")
#met = c("ald", "ald2", "ald5", "des", "des0.5", "des1")
#analysis.we[,3] <- rep(met,60)
    

# equivalence
# immuno.data_0.aldex2$p.data[[i]]$we.eBH == data.out[[i]]$ald0$wi.eBH
# immuno.data.DESeq$p.data[[i]]$padj == data.out[[i]]$des0$padj
#aldex.0 <- immuno.data_0.aldex2


row=1
for(coeff in c(0.01,0.1,0.2,0.5,0.75,1)){
for(i in 1:2){ 
# coefficients same in every instance
# can hardcode this if sanity check passes

model <- which(abs(input$thin.data[[11]]$coefmat) > coeff) 
null.model <- which(abs(input$thin.data[[11]]$coefmat) < coeff)


#for randomized data only (r)
#sample
TP.r <- intersect(which(input$r.data[[i]][,padj] < 0.05), model)
FN.r <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05))
FP.r <- intersect(which(input$r.data[[i]][,padj] < 0.05), null.model)
TN.r <- intersect(which(input$r.data[[i]][,padj] >= 0.05), null.model)

#for randomized + permuted data (p)
#sample
TP.p <- intersect(which(input$p.data[[i]][,padj] < 0.05), model)
FN.p <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05))
FP.p <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model)
TN.p <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model)

#how to integrate this
#TP.des5.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) >0.5), model)
#TP.des1.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) >1), model) # add in FC > 1
#FN.des5.r <- setdiff(model, which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 0.5))
#FN.des1.r <- setdiff(model, which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 1))
#FP.des5.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 0.5), null.model)
#FP.des1.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 1), null.model)
#TN.des5.r <- intersect(which(analysis.deseq$r.data[[i]]$padj >= 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 0.5), null.model)
#TN.des1.r <- intersect(which(analysis.deseq$r.data[[i]]$padj >= 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 1), null.model)

#for randomized data (r)
PPV.r <- length(TP.r)/sum(length(TP.r),length(FP.r))
FDR.r <- length(FP.r)/sum(length(TP.r),length(FP.r))
SEN.r <- length(TP.r)/(length(TP.r) + length(FN.r))
SPE.r <- length(TN.r)/(length(TN.r) + length(FP.r))

#for randomized + permuted data (p)
PPV.p <- length(TP.p)/sum(length(TP.p),length(FP.p))
FDR.p <- length(FP.p)/sum(length(TP.p),length(FP.p))
SEN.p <- length(TP.p)/(length(TP.p) + length(FN.p))
SPE.p <- length(TN.p)/(length(TN.p) + length(FP.p))

return(list(t_pos.r = TP.r, f_neg.r = FN.r, f_pos.r = FP.r, t_neg.r = TN.r))

} #end bracket here - finished editing up to here

# even butt uglier ...
for (j in 1:5){
  analysis.wi[row,1] <- coeff
  analysis.wi[row,2] <- i
if (j==1){
  analysis.wi[row,4] <- PPV.wi.ald
  analysis.wi[row,5] <- FDR.wi.ald
  analysis.wi[row,6] <- SEN.wi.ald
  analysis.wi[row,7] <- SPE.wi.ald
} else if (j==2) {
  analysis.wi[row,4] <- PPV.wi.ald.2
  analysis.wi[row,5] <- FDR.wi.ald.2
  analysis.wi[row,6] <- SEN.wi.ald.2
  analysis.wi[row,7] <- SPE.wi.ald.2
} else if (j==3){
  analysis.wi[row,4] <- PPV.wi.ald.5
  analysis.wi[row,5] <- FDR.wi.ald.5
  analysis.wi[row,6] <- SEN.wi.ald.5
  analysis.wi[row,7] <- SPE.wi.ald.5
} else if (j==4){
  analysis.wi[row,4] <- PPV.wi.des
  analysis.wi[row,5] <- FDR.wi.des
  analysis.wi[row,6] <- SEN.wi.des
  analysis.wi[row,7] <- SPE.wi.des
} else if (j==5){
  analysis.wi[row,4] <- PPV.wi.des5
  analysis.wi[row,5] <- FDR.wi.des5
  analysis.wi[row,6] <- SEN.wi.des5
  analysis.wi[row,7] <- SPE.wi.des5
}
row=row+1
}
}
}

means.wi <- matrix(data=NA, nrow=30, ncol=6)
means.wi <- as.data.frame(means.wi)
colnames(means.wi) <- c("coef", "method", "PPV","FDR","SEN","SPE")

# don't do this!!!!!!
# ugly hack
for(i in 1:5){
  means.wi[i,1:2] <- analysis.wi[seq(from=i, to=50, by=5),c(1,3)][1,]
  means.wi[i,3:6] <- colMeans(analysis.wi[seq(from=i, to=50, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+5,1:2] <- analysis.wi[seq(from=i+50, to=100, by=5),c(1,3)][1,]
  means.wi[i+5,3:6] <- colMeans(analysis.wi[seq(from=i+50, to=100, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+10,1:2] <- analysis.wi[seq(from=i+100, to=150, by=5),c(1,3)][1,]
  means.wi[i+10,3:6] <- colMeans(analysis.wi[seq(from=i+100, to=150, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+15,1:2] <- analysis.wi[seq(from=i+150, to=200, by=5),c(1,3)][1,]
  means.wi[i+15,3:6] <- colMeans(analysis.wi[seq(from=i+150, to=200, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+20,1:2] <- analysis.wi[seq(from=i+201, to=250, by=5),c(1,3)][1,]
  means.wi[i+20,3:6] <- colMeans(analysis.wi[seq(from=i+201, to=250, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+25,1:2] <- analysis.wi[seq(from=i+251, to=300, by=5),c(1,3)][1,]
  means.wi[i+25,3:6] <- colMeans(analysis.wi[seq(from=i+251, to=300, by=5),4:7])
}

FN.ald <- setdiff(model, which(data.out[[i]]$ald0$we.eBH < 0.05))
FN.ald.5 <- setdiff(model, which(data.out[[i]]$ald5$we.eBH < 0.05))
FN.ald.2 <- setdiff(model, which(data.out[[i]]$ald2$we.eBH < 0.05))
FN.des <- setdiff(model, which(data.out[[i]]$des$padj < 0.05))
FN.des5 <- setdiff(model, which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5))

FP.ald <- intersect(which(data.out[[i]]$ald0$we.eBH < 0.05), null.model)
FP.ald.2 <- intersect(which(data.out[[i]]$ald2$we.eBH < 0.05), null.model)
FP.ald.5 <- intersect(which(data.out[[i]]$ald5$we.eBH < 0.05), null.model)
FP.des <- intersect(which(data.out[[i]]$des$padj < 0.05), null.model)
FP.des5 <- intersect(which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)


TN.des <- intersect(which(data.out[[i]]$des$padj >= 0.05), null.model)
TN.des5 <- intersect(which(data.out[[i]]$des$padj >= 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)

PPV.ald <- length(TP.ald)/sum(length(TP.ald),length(FP.ald))
PPV.ald.2 <- length(TP.ald.2)/sum(length(TP.ald.2),length(FP.ald.2))
PPV.ald.5 <- length(TP.ald.5)/sum(length(TP.ald.5),length(FP.ald.5))
PPV.des <- length(TP.des)/sum(length(TP.des),length(FP.des))
PPV.des5 <- length(TP.des5)/sum(length(TP.des5),length(FP.des5))

FDR.ald <- length(FP.ald)/sum(length(TP.ald),length(FP.ald))
FDR.ald.2 <- length(FP.ald.2)/sum(length(TP.ald.2),length(FP.ald.2))
FDR.ald.5 <- length(FP.ald.5)/sum(length(TP.ald.5),length(FP.ald.5))
FDR.des <- length(FP.des)/sum(length(TP.des),length(FP.des))
FDR.des5 <- length(FP.des5)/sum(length(TP.des5),length(FP.des5))

SEN.ald <- length(TP.ald)/(length(TP.ald) + length(FN.ald))
SEN.ald.2 <- length(TP.ald.2)/(length(TP.ald.2) + length(FN.ald.2))
SEN.ald.5 <- length(TP.ald.5)/(length(TP.ald.5) + length(FN.ald.5))
SEN.des <- length(TP.des)/(length(TP.des) + length(FN.des))
SEN.des5 <- length(TP.des5)/(length(TP.des5) + length(FN.des5))

SPE.ald <- length(TN.ald)/(length(TN.ald) + length(FP.ald))
SPE.ald.2 <- length(TN.ald.2)/(length(TN.ald.2) + length(FP.ald.2))
SPE.ald.5 <- length(TN.ald.5)/(length(TN.ald.5) + length(FP.ald.5))
SPE.des <- length(TN.des)/(length(TN.des) + length(FP.des))
SPE.des5 <- length(TN.des5)/(length(TN.des5) + length(FP.des5))

# even butt uglier ...
for (j in 1:5){
  analysis.out[row,1] <- coeff
  analysis.out[row,2] <- i
if (j==1){
  analysis.out[row,4] <- PPV.ald
  analysis.out[row,5] <- FDR.ald
  analysis.out[row,6] <- SEN.ald
  analysis.out[row,7] <- SPE.ald
} else if (j==2) {
  analysis.out[row,4] <- PPV.ald.2
  analysis.out[row,5] <- FDR.ald.2
  analysis.out[row,6] <- SEN.ald.2
  analysis.out[row,7] <- SPE.ald.2
} else if (j==3){
  analysis.out[row,4] <- PPV.ald.5
  analysis.out[row,5] <- FDR.ald.5
  analysis.out[row,6] <- SEN.ald.5
  analysis.out[row,7] <- SPE.ald.5
} else if (j==4){
  analysis.out[row,4] <- PPV.des
  analysis.out[row,5] <- FDR.des
  analysis.out[row,6] <- SEN.des
  analysis.out[row,7] <- SPE.des
} else if (j==5){
  analysis.out[row,4] <- PPV.des5
  analysis.out[row,5] <- FDR.des5
  analysis.out[row,6] <- SEN.des5
  analysis.out[row,7] <- SPE.des5
}
row=row+1
}
}
}

means <- matrix(data=NA, nrow=30, ncol=6)
means <- as.data.frame(means)
colnames(means) <- c("coef", "method", "PPV","FDR","SEN","SPE")

for(i in 1:5){
  means[i,1:2] <- analysis.out[seq(from=i, to=50, by=5),c(1,3)][1,]
  means[i,3:6] <- colMeans(analysis.out[seq(from=i, to=50, by=5),4:7])
}
for(i in 1:5){
  means[i+5,1:2] <- analysis.out[seq(from=i+50, to=100, by=5),c(1,3)][1,]
  means[i+5,3:6] <- colMeans(analysis.out[seq(from=i+50, to=100, by=5),4:7])
}
for(i in 1:5){
  means[i+10,1:2] <- analysis.out[seq(from=i+100, to=150, by=5),c(1,3)][1,]
  means[i+10,3:6] <- colMeans(analysis.out[seq(from=i+100, to=150, by=5),4:7])
}
for(i in 1:5){
  means[i+15,1:2] <- analysis.out[seq(from=i+150, to=200, by=5),c(1,3)][1,]
  means[i+15,3:6] <- colMeans(analysis.out[seq(from=i+150, to=200, by=5),4:7])
}
for(i in 1:5){
  means[i+20,1:2] <- analysis.out[seq(from=i+201, to=250, by=5),c(1,3)][1,]
  means[i+20,3:6] <- colMeans(analysis.out[seq(from=i+201, to=250, by=5),4:7])
}
for(i in 1:5){
  means[i+25,1:2] <- analysis.out[seq(from=i+251, to=300, by=5),c(1,3)][1,]
  means[i+25,3:6] <- colMeans(analysis.out[seq(from=i+251, to=300, by=5),4:7])
}


ald.row <- which(means$method == "ald")
ald2.row <- which(means$method == "ald2")
ald5.row <- which(means$method == "ald5")
des.row <- which(means$method == "des")
des5.row <- which(means$method == "des5")
