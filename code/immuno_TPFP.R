analysis.fun <- function(input=NULL, type=NULL, nloop=100) {

  if (type == "DESeq") {
    padj = 6
    log2FoldChange = 3
  } else if (type == "edgeR") {
    padj = 5
    log2FoldChange = 2
  } else if (type == "aldex0.we") {
    padj = 10 #for welch's
    log2FoldChange = 5
  } else if (type == "aldex0.wi") {
    padj = 12  #for wilcoxon
    log2FoldChange = 5
  } else if (type == "aldex0.2.we") {
    padj = 10 
    log2FoldChange = 5
  } else if (type == "aldex0.2.wi") {
    padj = 12 
    log2FoldChange = 5
  } else if (type == "aldex0.5.we") {
    padj = 10 
    log2FoldChange = 5
  } else if (type == "aldex0.5.wi") {
    padj = 12
    log2FoldChange = 5
  } else if (type == "limma") {
    padj = 6 #use adj pvalue or just p value?
    log2FoldChange = 2
  }

#if (input == NULL) {stop("input modelled data")} #says to use is.null?
    
# equivalence
# immuno.data_0.aldex2$p.data[[i]]$we.eBH == data.out[[i]]$ald0$wi.eBH
# immuno.data.DESeq$p.data[[i]]$padj == data.out[[i]]$des0$padj
#aldex.0 <- immuno.data_0.aldex2

for(coeff in c(0.01,0.1,0.2,0.5,0.75,1)){

  for(i in 1:nloop){ 
  
  model <- list()
  null.model.r <- list()
  null.model.p <- list()
  TP.r <- list()
  FN.r <- list()
  FP.r <- list()
  TN.r <- list()
  
  TP.LFC.0.5.r <- list()
  TP.LFC.1.r <- list()
  FN.LFC.0.5.r <- list()
  FN.LFC.1.r <- list()
  FP.LFC.0.5.r <- list()
  FP.LFC.1.r <- list()
  TN.LFC.0.5.r <- list()
  TN.LFC.1.r <- list()
  
  TP.p <- list()
  FN.p <- list()
  FP.pr <- list()
  FP.pp <- list()
  TN.pr <- list()
  TN.pp <- list()
  
  TP.LFC.0.5.p <- list()
  TP.LFC.1.p <- list()
  FN.LFC.0.5.p <- list()
  FN.LFC.1.p <- list()
  FP.LFC.0.5.pr <- list()
  FP.LFC.1.pr <- list()
  FP.LFC.0.5.pp <- list()
  FP.LFC.1.pp <- list()
  TN.LFC.0.5.pr <- list()
  TN.LFC.1.pr <- list()
  TN.LFC.0.5.pp <- list()
  TN.LFC.1.pp <- list()
  
  PPV.r <- list()
  FDR.r <- list()
  SEN.r <- list()
  SPE.r <- list()
  
  PPV.pr <- list()
  PPV.pp <- list()
  FDR.pr <- list()
  FDR.pp <- list()
  SEN.p <- list()
  SPE.pr <- list()
  SPE.pp <- list()
  
# coefficients same in every instance
# can hardcode this if sanity check passes

model [[i]] <- which(abs(input$thin.data[[i]]$coefmat) > coeff)  
null.model.r [[i]] <- which(abs(input$thin.data[[i]]$coefmat) == 0) # randomized nulls only 
null.model.p [[i]] <- which(abs(input$thin.data[[i]]$coefmat) > 0 & abs(input$thin.data[[i]]$coefmat) < coeff) # only less than coeff FP 

#for randomized data only (r) 
TP.r[[i]] <- intersect(which(input$r.data[[i]][,padj] < 0.05), model)
FN.r[[i]] <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05))
FP.r[[i]] <- intersect(which(input$r.data[[i]][,padj] < 0.05), null.model.r) #uses only the randomized null model
TN.r[[i]] <- intersect(which(input$r.data[[i]][,padj] >= 0.05), null.model.r)

#for randomized data only (r) - DESeq with Log2FoldChange for 0.5 and 1, changed name to LFC, unsure if should rename?
TP.LFC.0.5.r[[i]] <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) >0.5), model)
TP.LFC.1.r[[i]] <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) >1), model) 
FN.LFC.0.5.r[[i]] <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5))
FN.LFC.1.r[[i]] <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1))
FP.LFC.0.5.r[[i]] <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
FP.LFC.1.r[[i]] <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1), null.model.r)
TN.LFC.0.5.r[[i]] <- intersect(which(input$r.data[[i]][,padj] >= 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
TN.LFC.1.r[[i]] <- intersect(which(input$r.data[[i]][,padj] >= 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1), null.model.r)

#for randomized + permuted data (p) 
# TN and FP uses both null models for randomized data and randomized+permuted data
TP.p[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05), model)
FN.p[[i]] <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05))
FP.pr[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model.r)
TN.pr[[i]] <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model.r)
FP.pp[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model.p)
TN.pp[[i]] <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model.p)

#for randomized data only (p) - DESeq with Log2FoldChange for 0.5 and 1
TP.LFC.0.5.p[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) >0.5), model)
TP.LFC.1.p[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) >1), model) 
FN.LFC.0.5.p[[i]] <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5))
FN.LFC.1.p[[i]] <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1))
FP.LFC.0.5.pr[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
FP.LFC.1.pr[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.r)
FP.LFC.0.5.pp[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.p)
FP.LFC.1.pp[[i]] <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.p)
TN.LFC.0.5.pr[[i]] <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
TN.LFC.1.pr[[i]] <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.r)
TN.LFC.0.5.pp[[i]] <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.p)
TN.LFC.1.pp[[i]] <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.p)

#for randomized data (r)
PPV.r[[i]] <- length(TP.r)/sum(length(TP.r),length(FP.r))
FDR.r[[i]] <- length(FP.r)/sum(length(TP.r),length(FP.r))
SEN.r[[i]] <- length(TP.r)/(length(TP.r) + length(FN.r))
SPE.r[[i]] <- length(TN.r)/(length(TN.r) + length(FP.r))

#for randomized + permuted data (p)
PPV.pr[[i]] <- length(TP.p)/sum(length(TP.p),length(FP.pr))
PPV.pp[[i]] <- length(TP.p)/sum(length(TP.p),length(FP.pp))
FDR.pr[[i]] <- length(FP.pr)/sum(length(TP.p),length(FP.pr))
FDR.pp[[i]] <- length(FP.pp)/sum(length(TP.p),length(FP.pp))
SEN.p[[i]] <- length(TP.p)/(length(TP.p) + length(FN.p)) #same for both?
SPE.pr[[i]] <- length(TN.pr)/(length(TN.pr) + length(FP.pr))
SPE.pp[[i]] <- length(TN.pp)/(length(TN.pp) + length(FP.pp))

return(list(truepos.r = TP.r, falseneg.r = FN.r, falsepos.r = FP.r, trueneg.r = TN.r, 
            TP.LFC.0.5.r = TP.LFC.0.5.r, TP.LFC.1.r = TP.LFC.1.r, FN.LFC.0.5.r = FN.LFC.0.5.r,  FN.LFC.1.r = FN.LFC.1.r, FP.LFC.0.5.r = FP.LFC.0.5.r, FP.LFC.1.r = FP.LFC.1.r, TN.LFC.0.5.r =  TN.LFC.0.5.r,  TN.LFC.1.r =  TN.LFC.1.r, 
            truepos.p = TP.p, falseneg.p = FN.r, falsepos.pr = FP.pr, falsepos.pp = FP.pp, trueneg.pr = TN.pr, trueneg.pp = TN.pp, 
            TP.LFC.0.5.p = TP.LFC.0.5.p, TP.LFC.1.p = TP.LFC.1.p, FN.LFC.0.5.p = FN.LFC.0.5.p,  FN.LFC.1.p = FN.LFC.1.p, FP.LFC.0.5.pr = FP.LFC.0.5.pr, FP.LFC.1.pr = FP.LFC.1.pr, TN.LFC.0.5.pr =  TN.LFC.0.5.pr,  TN.LFC.1.pr =  TN.LFC.1.pr, FP.LFC.0.5.pp = FP.LFC.0.5.pp, FP.LFC.1.pp = FP.LFC.1.pp, TN.LFC.0.5.pp =  TN.LFC.0.5.pp,  TN.LFC.1.pp =  TN.LFC.1.pp,
            PPV.r = PPV.r, FDR.r = FDR.r, SEN.r = SEN.r, SPE.r = SPE.r,
            PPV.pr = PPV.pr, PPV.pp = PPV.pp, FDR.pr = FDR.pr, FDR.pp = FDR.pp, SEN.p = SEN.p, SPE.pr = SPE.pr, SPE.pp = SPE.pp))
}}} #end bracket here - finished editing up to here

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
