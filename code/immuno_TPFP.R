load("immuno.data.edger.Rda")
load("immuno.data.deseq.Rda")
load("immuno.data.aldex2_0.out.Rda")
load("immuno.data.aldex2_2.out.Rda")
load("immuno.data.aldex2_5.out.Rda")

# "immuno.data_0.aldex2" "immuno.data_2.aldex2" "immuno.data_5.aldex2" "immuno.data.DESeq" "immuno.data.edgeR"  

analysis.fun <- function(input=NULL, type=NULL, nloop=100) {

  if (type = "DESeq") {
    padj = 6
  } else if (type = "edgeR") {
    padj = 5
  } else if (type = "aldex0.we") {
    we.eBH = 10
  } else if (type = "aldex0.wi") {
    wi.eBH = 12
  } else if (type = "aldex0.2.we"){
    we.eBH = 10
  } else if (type = "aldex0.2.wi"{
    wi.eBH = 12
  } else if (type = "aldex0.5.we") {
    we.eBH = 10
  } else if (type = "aldex0.5.wi") {
    wi.eBH = 12
  }
  
if(input == NULL){stop("input modelled data")}


### wilcoxon (t-test is labled welches below)
# there is duplication for deseq with the we and wilcoxon
analysis.wi <- matrix(data=NA, nrow=300, ncol=7) #change nrow/ncol?
analysis.wi <- as.data.frame(analysis.wi)
colnames(analysis.wi) <- c("coeff","iter", "met","PPV","FDR","SEN","SPE")
met = c("ald", "ald2", "ald5", "des", "des0.5", "des1")
analysis.wi[,3] <- rep(met,60)

#### we.eBH
analysis.we <- matrix(data=NA, nrow=300, ncol=7)
analysis.we <- as.data.frame(analysis.we)
colnames(analysis.we) <- c("coeff","iter", "met","PPV","FDR","SEN","SPE")
met = c("ald", "ald2", "ald5", "des", "des0.5", "des1")
analysis.we[,3] <- rep(met,60)
    

# equivalence
# immuno.data_0.aldex2$p.data[[i]]$we.eBH == data.out[[i]]$ald0$wi.eBH
# immuno.data.DESeq$p.data[[i]]$padj == data.out[[i]]$des0$padj
aldex.0 <- immuno.data_0.aldex2


row=1
for(coeff in c(0.01,0.1,0.2,0.5,0.75,1)){
for(i in 1:nloop){ 
# coefficients same in every instance
# can hardcode this if sanity check passes

model <- which(abs(analysis.edgeR$thin.data[[11]]$coefmat) > coeff) 
null.model <- which(abs(analysis.edgeR$thin.data[[11]]$coefmat) < coeff)


#for randomized data only (r)

#sample

TP.wi.ald.r  <- intersect(which(dataset$r.data[[11]]$wi.eBH < 0.05), model)
TP.wi.ald.5.r <- intersect(which(analysis.aldex0.5$r.data[[i]]$wi.eBH < 0.05), model)
TP.wi.ald.2.r <- intersect(which(analysis.aldex0.2$r.data[[i]]$wi.eBH < 0.05), model)
TP.we.ald.r  <- intersect(which(analysis.aldex0$r.data[[i]]$we.eBH < 0.05), model)
TP.we.ald.5.r <- intersect(which(analysis.aldex0.5$r.data[[i]]$we.eBH < 0.05), model)
TP.we.ald.2.r <- intersect(which(analysis.aldex0.2$r.data[[i]]$we.eBH < 0.05), model)
TP.wi.des.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05), model)
TP.des5.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) >0.5), model)
TP.des1.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) >1), model) # add in FC > 1

FN.wi.ald.r <- setdiff(model, which(analysis.aldex0$r.data[[i]]$wi.eBH < 0.05))
FN.wi.ald.5.r <- setdiff(model, which(analysis.aldex0.5$r.data[[i]]$wi.eBH < 0.05))
FN.wi.ald.2.r <- setdiff(model, which(analysis.aldex0.2$r.data[[i]]$wi.eBH < 0.05))
FN.we.ald.r <- setdiff(model, which(analysis.aldex0$r.data[[i]]$we.eBH < 0.05))
FN.we.ald.5.r <- setdiff(model, which(analysis.aldex0.5$r.data[[i]]$we.eBH < 0.05))
FN.we.ald.2.r <- setdiff(model, which(analysis.aldex0.2$r.data[[i]]$we.eBH < 0.05))
FN.des.r <- setdiff(model, which(analysis.deseq$r.data[[i]]$padj < 0.05))
FN.des5.r <- setdiff(model, which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 0.5))
FN.des1.r <- setdiff(model, which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 1))

FP.wi.ald.r <- intersect(which(analysis.aldex0$r.data[[i]]$wi.eBH < 0.05), null.model)
FP.wi.ald.2.r <- intersect(which(analysis.aldex0.2$r.data[[i]]$wi.eBH < 0.05), null.model)
FP.wi.ald.5.r <- intersect(which(analysis.aldex0.5$r.data[[i]]$wi.eBH < 0.05), null.model)
FP.we.ald.r <- intersect(which(analysis.aldex0$r.data[[i]]$we.eBH < 0.05), null.model)
TN.we.ald.2.r <- intersect(which(analysis.aldex0.2$r.data[[i]]$we.eBH < 0.05), null.model)
TN.we.ald.5.r <- intersect(which(analysis.aldex0.5$r.data[[i]]$we.eBH < 0.05), null.model)
FP.des.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05), null.model)
FP.des5.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 0.5), null.model)
FP.des1.r <- intersect(which(analysis.deseq$r.data[[i]]$padj < 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 1), null.model)

TN.wi.ald.r <- intersect(which(analysis.aldex0$r.data[[i]]$wi.eBH >= 0.05), null.model)
TN.wi.ald.2.r <- intersect(which(analysis.aldex0.2$r.data[[i]]$wi.eBH >= 0.05), null.model)
TN.wi.ald.5.r <- intersect(which(analysis.aldex0.5$r.data[[i]]$wi.eBH >= 0.05), null.model)
TN.we.ald.r <- intersect(which(analysis.aldex0$r.data[[i]]$we.eBH >= 0.05), null.model)
TN.we.ald.2.r <- intersect(which(analysis.aldex0.2$r.data[[i]]$we.eBH >= 0.05), null.model)
TN.we.ald.5.r <- intersect(which(analysis.aldex0.5$r.data[[i]]$we.eBH >= 0.05), null.model)
TN.des.r <- intersect(which(analysis.deseq$r.data[[i]]$padj >= 0.05), null.model)
TN.des5.r <- intersect(which(analysis.deseq$r.data[[i]]$padj >= 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 0.5), null.model)
TN.des1.r <- intersect(which(analysis.deseq$r.data[[i]]$padj >= 0.05 & abs(analysis.deseq$r.data[[i]]$log2FoldChange) > 1), null.model)


#for permuted + randomized data (p)
TP.wi.ald.p  <- intersect(which(analysis.aldex0$p.data[[11]]$wi.eBH < 0.05), model)
TP.wi.ald.5.p <- intersect(which(analysis.aldex0.5$p.data[[i]]$wi.eBH < 0.05), model)
TP.wi.ald.2.p <- intersect(which(analysis.aldex0.2$p.data[[i]]$wi.eBH < 0.05), model)
TP.we.ald.p  <- intersect(which(analysis.aldex0$p.data[[i]]$we.eBH < 0.05), model)
TP.we.ald.5.p <- intersect(which(analysis.aldex0.5$p.data[[i]]$we.eBH < 0.05), model)
TP.we.ald.2.p <- intersect(which(analysis.aldex0.2$p.data[[i]]$we.eBH < 0.05), model)
TP.des.p <- intersect(which(analysis.deseq$p.data[[i]]$padj < 0.05), model)
TP.des5.p <- intersect(which(analysis.deseq$p.data[[i]]$padj < 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) >0.5), model)
TP.des1.p <- intersect(which(analysis.deseq$p.data[[i]]$padj < 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) >1), model) # add in FC > 1

FN.wi.ald.p <- setdiff(model, which(analysis.aldex0$p.data[[i]]$wi.eBH < 0.05))
FN.wi.ald.5.p <- setdiff(model, which(analysis.aldex0.5$p.data[[i]]$wi.eBH < 0.05))
FN.wi.ald.2.p <- setdiff(model, which(analysis.aldex0.2$p.data[[i]]$wi.eBH < 0.05))
FN.we.ald.p <- setdiff(model, which(analysis.aldex0$p.data[[i]]$we.eBH < 0.05))
FN.we.ald.5.p <- setdiff(model, which(analysis.aldex0.5$p.data[[i]]$we.eBH < 0.05))
FN.we.ald.2.p <- setdiff(model, which(analysis.aldex0.2$p.data[[i]]$we.eBH < 0.05))
FN.des.p <- setdiff(model, which(analysis.deseq$p.data[[i]]$padj < 0.05))
FN.des5.p <- setdiff(model, which(analysis.deseq$p.data[[i]]$padj < 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) > 0.5))
FN.des1.p <- setdiff(model, which(analysis.deseq$p.data[[i]]$padj < 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) > 1))

FP.wi.ald.p <- intersect(which(analysis.aldex0$p.data[[i]]$wi.eBH < 0.05), null.model)
FP.wi.ald.2.p <- intersect(which(analysis.aldex0.2$p.data[[i]]$wi.eBH < 0.05), null.model)
FP.wi.ald.5.p <- intersect(which(analysis.aldex0.5$p.data[[i]]$wi.eBH < 0.05), null.model)
FP.we.ald.p <- intersect(which(analysis.aldex0$p.data[[i]]$we.eBH < 0.05), null.model)
TN.we.ald.2.p <- intersect(which(analysis.aldex0.2$p.data[[i]]$we.eBH < 0.05), null.model)
TN.we.ald.5.p <- intersect(which(analysis.aldex0.5$p.data[[i]]$we.eBH < 0.05), null.model)
FP.des.p <- intersect(which(analysis.deseq$p.data[[i]]$padj < 0.05), null.model)
FP.des5.p <- intersect(which(analysis.deseq$p.data[[i]]$padj < 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) > 0.5), null.model)
FP.des1.p <- intersect(which(analysis.deseq$p.data[[i]]$padj < 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) > 1), null.model)

TN.wi.ald.p <- intersect(which(analysis.aldex0$p.data[[i]]$wi.eBH >= 0.05), null.model)
TN.wi.ald.2.p <- intersect(which(analysis.aldex0.2$p.data[[i]]$wi.eBH >= 0.05), null.model)
TN.wi.ald.5.p <- intersect(which(analysis.aldex0.5$p.data[[i]]$wi.eBH >= 0.05), null.model)
TN.we.ald.p <- intersect(which(analysis.aldex0$p.data[[i]]$we.eBH >= 0.05), null.model)
TN.we.ald.2.p <- intersect(which(analysis.aldex0.2$p.data[[i]]$we.eBH >= 0.05), null.model)
TN.we.ald.5.p <- intersect(which(analysis.aldex0.5$p.data[[i]]$we.eBH >= 0.05), null.model)
TN.des.p <- intersect(which(analysis.deseq$p.data[[i]]$padj >= 0.05), null.model)
TN.des5.p <- intersect(which(analysis.deseq$p.data[[i]]$padj >= 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) > 0.5), null.model)
TN.des1.p <- intersect(which(analysis.deseq$p.data[[i]]$padj >= 0.05 & abs(analysis.deseq$p.data[[i]]$log2FoldChange) > 1), null.model)

# should not have to change but can augment with edgeR, limma, etc
#for randomized data (r)
PPV.wi.ald.r <- length(TP.wi.ald.r)/sum(length(TP.wi.ald.r),length(FP.wi.ald.r))
PPV.wi.ald.2.r <- length(TP.wi.ald.2.r)/sum(length(TP.wi.ald.2.r),length(FP.wi.ald.2.r))
PPV.wi.ald.5.r <- length(TP.wi.ald.5.r)/sum(length(TP.wi.ald.5.r),length(FP.wi.ald.5.r))
PPV.we.ald.r <- length(TP.we.ald.r)/sum(length(TP.we.ald.r),length(FP.we.ald.r))
PPV.we.ald.2.r <- length(TP.we.ald.2.r)/sum(length(TP.we.ald.2.r),length(FP.we.ald.2.r))
PPV.we.ald.5.r <- length(TP.we.ald.5.r)/sum(length(TP.we.ald.5.r),length(FP.we.ald.5.r))
PPV.des.r <- length(TP.des.r)/sum(length(TP.des.r),length(FP.des.r))
PPV.des5.r <- length(TP.des5.r)/sum(length(TP.des5.r),length(FP.des5.r))
PPV.des1.r <- length(TP.des1.r)/sum(length(TP.des1.r),length(FP.des1.r))


FDR.wi.ald.r <- length(FP.wi.ald.r)/sum(length(TP.wi.ald.r),length(FP.wi.ald.r))
FDR.wi.ald.2.r <- length(FP.wi.ald.2.r)/sum(length(TP.wi.ald.2.r),length(FP.wi.ald.2.r))
FDR.wi.ald.5.r <- length(FP.wi.ald.5.r)/sum(length(TP.wi.ald.5.r),length(FP.wi.ald.5.r))
FDR.we.ald.r <- length(FP.we.ald.r)/sum(length(TP.we.ald.r),length(FP.we.ald.r))
FDR.we.ald.2.r <- length(FP.we.ald.2.r)/sum(length(TP.we.ald.2.r),length(FP.we.ald.2.r))
FDR.we.ald.5.r <- length(FP.we.ald.5.r)/sum(length(TP.we.ald.5.r),length(FP.we.ald.5.r))
FDR.wi.des <- length(FP.wi.des.r)/sum(length(TP.wi.des.r),length(FP.wi.des.r))
FDR.wi.des5 <- length(FP.wi.des5.r)/sum(length(TP.wi.des5.r),length(FP.wi.des5.r))

SEN.wi.ald <- length(TP.wi.ald)/(length(TP.wi.ald) + length(FN.wi.ald))
SEN.wi.ald.2 <- length(TP.wi.ald.2)/(length(TP.wi.ald.2) + length(FN.wi.ald.2))
SEN.wi.ald.5 <- length(TP.wi.ald.5)/(length(TP.wi.ald.5) + length(FN.wi.ald.5))
SEN.wi.des <- length(TP.wi.des)/(length(TP.wi.des) + length(FN.wi.des))
SEN.wi.des5 <- length(TP.wi.des5)/(length(TP.wi.des5) + length(FN.wi.des5))

SPE.wi.ald <- length(TN.wi.ald)/(length(TN.wi.ald) + length(FP.wi.ald))
SPE.wi.ald.2 <- length(TN.wi.ald.2)/(length(TN.wi.ald.2) + length(FP.wi.ald.2))
SPE.wi.ald.5 <- length(TN.wi.ald.5)/(length(TN.wi.ald.5) + length(FP.wi.ald.5))
SPE.wi.des <- length(TN.wi.des)/(length(TN.wi.des) + length(FP.wi.des))
SPE.wi.des5 <- length(TN.wi.des5)/(length(TN.wi.des5) + length(FP.wi.des5))

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
