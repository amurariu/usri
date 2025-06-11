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

	coeff.vec <- c(0.01,0.1,0.2,0.5,0.75,1)
		
	PPV.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SEN.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	PPV.pr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	PPV.pp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.pr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.pp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SEN.p <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.pr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.pp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	for(coeff in 1:6){
			for(i in 1:nloop){ 
			
			
		# coefficients same in every instance
		# can hardcode this if sanity check passes
		
		model <- which(abs(input$thin.data[[i]]$coefmat) > coeff.vec[coeff])  
		null.model.r <- which(abs(input$thin.data[[i]]$coefmat) == 0) # randomized nulls only 
		null.model.p <- which(abs(input$thin.data[[i]]$coefmat) > 0 & abs(input$thin.data[[i]]$coefmat) < coeff) # only less than coeff FP 
		
		#for randomized data only (r) 
		TP.r <- intersect(which(input$r.data[[i]][,padj] < 0.05), model)
		FN.r <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05))
		FP.r <- intersect(which(input$r.data[[i]][,padj] < 0.05), null.model.r) #uses only the randomized null model
		TN.r <- intersect(which(input$r.data[[i]][,padj] >= 0.05), null.model.r)
		
		#for randomized data only (r) - DESeq with Log2FoldChange for 0.5 and 1, changed name to LFC, unsure if should rename?
		TP.LFC.0.5.r <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) >0.5), model)
		TP.LFC.1.r <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) >1), model) 
		FN.LFC.0.5.r <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5))
		FN.LFC.1.r <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1))
		FP.LFC.0.5.r <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		FP.LFC.1.r <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1), null.model.r)
		TN.LFC.0.5.r <- intersect(which(input$r.data[[i]][,padj] >= 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		TN.LFC.1.r <- intersect(which(input$r.data[[i]][,padj] >= 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1), null.model.r)
		
		#for randomized + permuted data (p) 
		# TN and FP uses both null models for randomized data and randomized+permuted data
		TP.p <- intersect(which(input$p.data[[i]][,padj] < 0.05), model)
		FN.p <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05))
		FP.pr <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model.r)
		TN.pr <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model.r)
		FP.pp <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model.p)
		TN.pp <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model.p)
		
		#for randomized data only (p) - DESeq with Log2FoldChange for 0.5 and 1
		TP.LFC.0.5.p <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) >0.5), model)
		TP.LFC.1.p <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) >1), model) 
		FN.LFC.0.5.p <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5))
		FN.LFC.1.p <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1))
		FP.LFC.0.5.pr <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		FP.LFC.1.pr <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.r)
		FP.LFC.0.5.pp <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.p)
		FP.LFC.1.pp <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.p)
		TN.LFC.0.5.pr <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		TN.LFC.1.pr <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.r)
		TN.LFC.0.5.pp <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.p)
		TN.LFC.1.pp <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.p)
		
		#for randomized data (r)
		# if PPV.r is a matrix PPV[coeff,i]
		PPV.r[i, coeff] <- length(TP.r)/sum(length(TP.r),length(FP.r))
		FDR.r[i, coeff] <- length(FP.r)/sum(length(TP.r),length(FP.r))
		SEN.r[i, coeff] <- length(TP.r)/(length(TP.r) + length(FN.r))
		SPE.r[i, coeff] <- length(TN.r)/(length(TN.r) + length(FP.r))
		
		#for randomized + permuted data (p)
		PPV.pr[i, coeff] <- length(TP.p)/sum(length(TP.p),length(FP.pr))
		PPV.pp[i, coeff] <- length(TP.p)/sum(length(TP.p),length(FP.pp))
		FDR.pr[i, coeff] <- length(FP.pr)/sum(length(TP.p),length(FP.pr))
		FDR.pp[i, coeff] <- length(FP.pp)/sum(length(TP.p),length(FP.pp))
		SEN.p[i, coeff] <- length(TP.p)/(length(TP.p) + length(FN.p)) #same for both?
		SPE.pr[i, coeff] <- length(TN.pr)/(length(TN.pr) + length(FP.pr))
		SPE.pp[i, coeff] <- length(TN.pp)/(length(TN.pp) + length(FP.pp))
		
		return(list(PPV.r, FDR.r, SEN.r, SPE.r, PPV.pr,  PPV.pp, FDR.pr, FDR.pp, SEN.p, SPE.pr = SPE.pr, SPE.pp = SPE.pp))
		} # end coeff loop
	} # end randomization loop
} #end function 
# finished editing up to here

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
