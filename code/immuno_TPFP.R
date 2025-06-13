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

	coeff.vec <- c(0.01,0.1,0.2,0.5,0.75,1)
		
	PPV.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SEN.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	PPV.perm.nullr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	PPV.perm.nullp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.perm.nullr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.perm.nullp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SEN.perm <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.perm.nullr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.perm.nullp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	PPV.LFC.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.LFC.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SEN.LFC.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.LFC.rand <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	#PPV.LFC.1.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	#FDR.LFC.1.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	#SEN.LFC.1.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	#SPE.LFC.1.r <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	PPV.LFC.perm.nullp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	PPV.LFC.perm.nullr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	FDR.LFC.perm.nullp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	FDR.LFC.perm.nullr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	SEN.LFC.perm <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	SPE.LFC.perm.nullp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	SPE.LFC.perm.nullr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	#PPV.LFC.1.pp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	#PPV.LFC.1.pr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	#FDR.LFC.1.pp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	#FDR.LFC.1.pr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	#SEN.LFC.1.p <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	#SPE.LFC.1.pp <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	#SPE.LFC.1.pr <- matrix(data=NA, nrow=nloop, ncol=length(coeff.vec))
	
	for(coeff in 1:6){
			for(i in 1:nloop){ 
		 # print(c(coeff, i))
		# coefficients same in every instance
		# can hardcode this if sanity check passes
		
		model <- which(abs(input$thin.data[[i]]$coefmat) > coeff.vec[coeff])  
		null.model.r <- which(abs(input$thin.data[[i]]$coefmat) == 0) # randomized nulls only 
		null.model.p <- which(abs(input$thin.data[[i]]$coefmat) > 0 & abs(input$thin.data[[i]]$coefmat) < coeff) # only less than coeff FP 
		
		#for randomized data only (r) 
		TP.rand <- intersect(which(input$r.data[[i]][,padj] < 0.05), model)
		FN.rand <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05))
		FP.rand <- intersect(which(input$r.data[[i]][,padj] < 0.05), null.model.r) #uses only the randomized null model
		TN.rand <- intersect(which(input$r.data[[i]][,padj] >= 0.05), null.model.r)
		
		#for randomized data only (r) - DESeq with Log2FoldChange for 0.5 and 1, changed name to LFC, unsure if should rename?
		TP.LFC.rand <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) >0.5), model)
		#TP.LFC.1.r <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) >1), model) #no entries/IDs identified
		FN.LFC.rand <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5))
		#FN.LFC.1.r <- setdiff(model, which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1))
		FP.LFC.rand <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		#FP.LFC.1.r <- intersect(which(input$r.data[[i]][,padj] < 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1), null.model.r) #no entries/IDs identified
		TN.LFC.rand <- intersect(which(input$r.data[[i]][,padj] >= 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		#TN.LFC.1.r <- intersect(which(input$r.data[[i]][,padj] >= 0.05 & abs(input$r.data[[i]][,log2FoldChange]) > 1), null.model.r)
		
		#for randomized + permuted data (p) 
		# TN and FP uses both null models for randomized data and randomized+permuted data
		TP.perm <- intersect(which(input$p.data[[i]][,padj] < 0.05), model)
		FN.perm <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05))
		FP.perm <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model.r)
		TN.perm <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model.r)
		FP.perm <- intersect(which(input$p.data[[i]][,padj] < 0.05), null.model.p)
		TN.perm <- intersect(which(input$p.data[[i]][,padj] >= 0.05), null.model.p)
		
		#for randomized data only (p) - DESeq with Log2FoldChange for 0.5 and 1
		TP.LFC.perm <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) >0.5), model)
		#TP.LFC.1.p <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) >1), model) #no entries/IDs identified
		FN.LFC.perm <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5))
		#FN.LFC.1.p <- setdiff(model, which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1))
		FP.LFC.perm.nullr <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		#FP.LFC.1.pr <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.r) #no entries/IDs identified
		FP.LFC.perm.nullp <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.p)
		#FP.LFC.1.pp <- intersect(which(input$p.data[[i]][,padj] < 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.p) #no entries/IDs identified
		TN.LFC.perm.nullr <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.r)
		#TN.LFC.1.pr <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.r)
		TN.LFC.perm.nullp <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 0.5), null.model.p)
		#TN.LFC.1.pp <- intersect(which(input$p.data[[i]][,padj] >= 0.05 & abs(input$p.data[[i]][,log2FoldChange]) > 1), null.model.p) #no entries/IDs identified
		
		#for randomized data (r)
		# if PPV.r is a matrix PPV[coeff,i]
		PPV.rand[i, coeff] <- length(TP.r)/sum(length(TP.r),length(FP.r))
		FDR.rand[i, coeff] <- length(FP.r)/sum(length(TP.r),length(FP.r))
		SEN.rand[i, coeff] <- length(TP.r)/(length(TP.r) + length(FN.r))
		SPE.rand[i, coeff] <- length(TN.r)/(length(TN.r) + length(FP.r))
		
		PPV.LFC.rand[i, coeff] <- length(TP.LFC.0.5.r)/sum(length(TP.LFC.0.5.r),length(FP.LFC.0.5.r))
		FDR.LFC.rand[i, coeff] <- length(FP.LFC.0.5.r)/sum(length(TP.LFC.0.5.r),length(FP.LFC.0.5.r))
		SEN.LFC.rand[i, coeff] <- length(TP.LFC.0.5.r)/(length(TP.LFC.0.5.r) + length(FN.LFC.0.5.r))
		SPE.LFC.rand[i, coeff] <- length(TN.LFC.0.5.r)/(length(TN.LFC.0.5.r) + length(FP.LFC.0.5.r))
		
		#PPV.LFC.1.r[i, coeff] <- length(TP.LFC.1.r)/sum(length(TP.LFC.1.r),length(FP.LFC.1.r)) #NaN
		#FDR.LFC.1.r[i, coeff] <- length(FP.LFC.1.r)/sum(length(TP.LFC.1.r),length(FP.LFC.1.r)) #NaN
		#SEN.LFC.1.r[i, coeff] <- length(TP.LFC.1.r)/(length(TP.LFC.1.r) + length(FN.LFC.1.r)) #0
		#SPE.LFC.1.r[i, coeff] <- length(TN.LFC.1.r)/(length(TN.LFC.1.r) + length(FP.LFC.1.r)) #1
		
		
		#for randomized + permuted data (p)
		PPV.perm.nullr[i, coeff] <- length(TP.p)/sum(length(TP.p),length(FP.pr))
		PPV.perm.nullp[i, coeff] <- length(TP.p)/sum(length(TP.p),length(FP.pp))
		FDR.perm.nullr[i, coeff] <- length(FP.pr)/sum(length(TP.p),length(FP.pr))
		FDR.perm.nullp[i, coeff] <- length(FP.pp)/sum(length(TP.p),length(FP.pp))
		SEN.perm[i, coeff] <- length(TP.p)/(length(TP.p) + length(FN.p)) 
		SPE.perm.nullr[i, coeff] <- length(TN.pr)/(length(TN.pr) + length(FP.pr))
		SPE.perm.nullp[i, coeff] <- length(TN.pp)/(length(TN.pp) + length(FP.pp))
		
		PPV.LFC.perm.nullp[i, coeff] <- length(TP.LFC.0.5.p)/sum(length(TP.LFC.0.5.p),length(FP.LFC.0.5.pp))
		PPV.LFC.perm.nullr[i, coeff] <- length(TP.LFC.0.5.p)/sum(length(TP.LFC.0.5.p),length(FP.LFC.0.5.pr))
		
		FDR.LFC.perm.nullp[i, coeff] <- length(FP.LFC.0.5.pp)/sum(length(TP.LFC.0.5.p),length(FP.LFC.0.5.pp))
		FDR.LFC.perm.nullr[i, coeff] <- length(FP.LFC.0.5.pr)/sum(length(TP.LFC.0.5.p),length(FP.LFC.0.5.pr))
		
		SEN.LFC.perm[i, coeff] <- length(TP.LFC.0.5.p)/(length(TP.LFC.0.5.p) + length(FN.LFC.0.5.p))
		
		SPE.LFC.perm.nullp[i, coeff] <- length(TN.LFC.0.5.pp)/(length(TN.LFC.0.5.pp) + length(FP.LFC.0.5.pp))
		SPE.LFC.perm.nullr[i, coeff] <- length(TN.LFC.0.5.pr)/(length(TN.LFC.0.5.pr) + length(FP.LFC.0.5.pr))
		
		
		#PPV.LFC.1.pp[i, coeff] <- length(TP.LFC.1.p)/sum(length(TP.LFC.1.p),length(FP.LFC.1.pp)) #NaN
		#PPV.LFC.1.pr[i, coeff] <- length(TP.LFC.1.p)/sum(length(TP.LFC.1.p),length(FP.LFC.1.pr)) #NaN
		
		#FDR.LFC.1.pp[i, coeff] <- length(FP.LFC.1.pp)/sum(length(TP.LFC.1.p),length(FP.LFC.1.pp)) #NaN
		#FDR.LFC.1.pr[i, coeff] <- length(FP.LFC.1.pr)/sum(length(TP.LFC.1.p),length(FP.LFC.1.pr)) #NaN
		
		#SEN.LFC.1.p[i, coeff] <- length(TP.LFC.1.p)/(length(TP.LFC.1.p) + length(FN.LFC.1.p)) #0
		
		#SPE.LFC.1.pp[i, coeff] <- length(TN.LFC.1.pp)/(length(TN.LFC.1.pp) + length(FP.LFC.1.pp)) #NaN
		#SPE.LFC.1.pr[i, coeff] <- length(TN.LFC.1.pr)/(length(TN.LFC.1.pr) + length(FP.LFC.1.pr)) #1


		} # end coeff loop
	} # end randomization loop
	return(list(TP.rand = TP.rand, FN.rand = FN.rand, FP.rand = FP.rand, TN.rand = TN.rand, PPV.rand = PPV.rand, FDR.rand = FDR.rand, SEN.rand = SEN.rand, SPE.rand = SPE.rand, PPV.perm.nullr = PPV.perm.nullr,  PPV.perm.nullp = PPV.perm.nullp, FDR.perm.nullr = FDR.perm.nullr, FDR.perm.nullp = FDR.perm.nullp, SEN.perm = SEN.perm, SPE.perm.nullr = SPE.perm.nullr, SPE.perm.nullp = SPE.perm.nullp, PPV.LFC.rand = PPV.LFC.rand, FDR.LFC.rand = FDR.LFC.rand, SPE.LFC.rand = SPE.LFC.rand, SEN.LFC.rand = SEN.LFC.rand, PPV.LFC.perm.nullp = PPV.LFC.perm.nullp, PPV.LFC.perm.nullr = PPV.LFC.perm.nullr, FDR.LFC.perm.nullp = FDR.LFC.perm.nullp, FDR.LFC.perm.nullr, SEN.LFC.perm = SEN.LFC.perm, SPE.LFC.perm.nullp = SPE.LFC.perm.nullp, SPE.LFC.perm.nullr = SPE.LFC.perm.nullr))
} #end function 
# finished editing up to here
# 
# # even butt uglier ...
# for (j in 1:5){
#   analysis.wi[row,1] <- coeff
#   analysis.wi[row,2] <- i
# if (j==1){
#   analysis.wi[row,4] <- PPV.wi.ald
#   analysis.wi[row,5] <- FDR.wi.ald
#   analysis.wi[row,6] <- SEN.wi.ald
#   analysis.wi[row,7] <- SPE.wi.ald
# } else if (j==2) {
#   analysis.wi[row,4] <- PPV.wi.ald.2
#   analysis.wi[row,5] <- FDR.wi.ald.2
#   analysis.wi[row,6] <- SEN.wi.ald.2
#   analysis.wi[row,7] <- SPE.wi.ald.2
# } else if (j==3){
#   analysis.wi[row,4] <- PPV.wi.ald.5
#   analysis.wi[row,5] <- FDR.wi.ald.5
#   analysis.wi[row,6] <- SEN.wi.ald.5
#   analysis.wi[row,7] <- SPE.wi.ald.5
# } else if (j==4){
#   analysis.wi[row,4] <- PPV.wi.des
#   analysis.wi[row,5] <- FDR.wi.des
#   analysis.wi[row,6] <- SEN.wi.des
#   analysis.wi[row,7] <- SPE.wi.des
# } else if (j==5){
#   analysis.wi[row,4] <- PPV.wi.des5
#   analysis.wi[row,5] <- FDR.wi.des5
#   analysis.wi[row,6] <- SEN.wi.des5
#   analysis.wi[row,7] <- SPE.wi.des5
# }
# row=row+1
# }
# }
# }
# 
# means.wi <- matrix(data=NA, nrow=30, ncol=6)
# means.wi <- as.data.frame(means.wi)
# colnames(means.wi) <- c("coef", "method", "PPV","FDR","SEN","SPE")
# 
# # don't do this!!!!!!
# # ugly hack
# for(i in 1:5){
#   means.wi[i,1:2] <- analysis.wi[seq(from=i, to=50, by=5),c(1,3)][1,]
#   means.wi[i,3:6] <- colMeans(analysis.wi[seq(from=i, to=50, by=5),4:7])
# }
# for(i in 1:5){
#   means.wi[i+5,1:2] <- analysis.wi[seq(from=i+50, to=100, by=5),c(1,3)][1,]
#   means.wi[i+5,3:6] <- colMeans(analysis.wi[seq(from=i+50, to=100, by=5),4:7])
# }
# for(i in 1:5){
#   means.wi[i+10,1:2] <- analysis.wi[seq(from=i+100, to=150, by=5),c(1,3)][1,]
#   means.wi[i+10,3:6] <- colMeans(analysis.wi[seq(from=i+100, to=150, by=5),4:7])
# }
# for(i in 1:5){
#   means.wi[i+15,1:2] <- analysis.wi[seq(from=i+150, to=200, by=5),c(1,3)][1,]
#   means.wi[i+15,3:6] <- colMeans(analysis.wi[seq(from=i+150, to=200, by=5),4:7])
# }
# for(i in 1:5){
#   means.wi[i+20,1:2] <- analysis.wi[seq(from=i+201, to=250, by=5),c(1,3)][1,]
#   means.wi[i+20,3:6] <- colMeans(analysis.wi[seq(from=i+201, to=250, by=5),4:7])
# }
# for(i in 1:5){
#   means.wi[i+25,1:2] <- analysis.wi[seq(from=i+251, to=300, by=5),c(1,3)][1,]
#   means.wi[i+25,3:6] <- colMeans(analysis.wi[seq(from=i+251, to=300, by=5),4:7])
# }
# 
# FN.ald <- setdiff(model, which(data.out[[i]]$ald0$we.eBH < 0.05))
# FN.ald.5 <- setdiff(model, which(data.out[[i]]$ald5$we.eBH < 0.05))
# FN.ald.2 <- setdiff(model, which(data.out[[i]]$ald2$we.eBH < 0.05))
# FN.des <- setdiff(model, which(data.out[[i]]$des$padj < 0.05))
# FN.des5 <- setdiff(model, which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5))
# 
# FP.ald <- intersect(which(data.out[[i]]$ald0$we.eBH < 0.05), null.model)
# FP.ald.2 <- intersect(which(data.out[[i]]$ald2$we.eBH < 0.05), null.model)
# FP.ald.5 <- intersect(which(data.out[[i]]$ald5$we.eBH < 0.05), null.model)
# FP.des <- intersect(which(data.out[[i]]$des$padj < 0.05), null.model)
# FP.des5 <- intersect(which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)
# 
# 
# TN.des <- intersect(which(data.out[[i]]$des$padj >= 0.05), null.model)
# TN.des5 <- intersect(which(data.out[[i]]$des$padj >= 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)
# 
# PPV.ald <- length(TP.ald)/sum(length(TP.ald),length(FP.ald))
# PPV.ald.2 <- length(TP.ald.2)/sum(length(TP.ald.2),length(FP.ald.2))
# PPV.ald.5 <- length(TP.ald.5)/sum(length(TP.ald.5),length(FP.ald.5))
# PPV.des <- length(TP.des)/sum(length(TP.des),length(FP.des))
# PPV.des5 <- length(TP.des5)/sum(length(TP.des5),length(FP.des5))
# 
# FDR.ald <- length(FP.ald)/sum(length(TP.ald),length(FP.ald))
# FDR.ald.2 <- length(FP.ald.2)/sum(length(TP.ald.2),length(FP.ald.2))
# FDR.ald.5 <- length(FP.ald.5)/sum(length(TP.ald.5),length(FP.ald.5))
# FDR.des <- length(FP.des)/sum(length(TP.des),length(FP.des))
# FDR.des5 <- length(FP.des5)/sum(length(TP.des5),length(FP.des5))
# 
# SEN.ald <- length(TP.ald)/(length(TP.ald) + length(FN.ald))
# SEN.ald.2 <- length(TP.ald.2)/(length(TP.ald.2) + length(FN.ald.2))
# SEN.ald.5 <- length(TP.ald.5)/(length(TP.ald.5) + length(FN.ald.5))
# SEN.des <- length(TP.des)/(length(TP.des) + length(FN.des))
# SEN.des5 <- length(TP.des5)/(length(TP.des5) + length(FN.des5))
# 
# SPE.ald <- length(TN.ald)/(length(TN.ald) + length(FP.ald))
# SPE.ald.2 <- length(TN.ald.2)/(length(TN.ald.2) + length(FP.ald.2))
# SPE.ald.5 <- length(TN.ald.5)/(length(TN.ald.5) + length(FP.ald.5))
# SPE.des <- length(TN.des)/(length(TN.des) + length(FP.des))
# SPE.des5 <- length(TN.des5)/(length(TN.des5) + length(FP.des5))
# 
# # even butt uglier ...
# for (j in 1:5){
#   analysis.out[row,1] <- coeff
#   analysis.out[row,2] <- i
# if (j==1){
#   analysis.out[row,4] <- PPV.ald
#   analysis.out[row,5] <- FDR.ald
#   analysis.out[row,6] <- SEN.ald
#   analysis.out[row,7] <- SPE.ald
# } else if (j==2) {
#   analysis.out[row,4] <- PPV.ald.2
#   analysis.out[row,5] <- FDR.ald.2
#   analysis.out[row,6] <- SEN.ald.2
#   analysis.out[row,7] <- SPE.ald.2
# } else if (j==3){
#   analysis.out[row,4] <- PPV.ald.5
#   analysis.out[row,5] <- FDR.ald.5
#   analysis.out[row,6] <- SEN.ald.5
#   analysis.out[row,7] <- SPE.ald.5
# } else if (j==4){
#   analysis.out[row,4] <- PPV.des
#   analysis.out[row,5] <- FDR.des
#   analysis.out[row,6] <- SEN.des
#   analysis.out[row,7] <- SPE.des
# } else if (j==5){
#   analysis.out[row,4] <- PPV.des5
#   analysis.out[row,5] <- FDR.des5
#   analysis.out[row,6] <- SEN.des5
#   analysis.out[row,7] <- SPE.des5
# }
# row=row+1
# }
# }
# }
# 
# means <- matrix(data=NA, nrow=30, ncol=6)
# means <- as.data.frame(means)
# colnames(means) <- c("coef", "method", "PPV","FDR","SEN","SPE")
# 
# for(i in 1:5){
#   means[i,1:2] <- analysis.out[seq(from=i, to=50, by=5),c(1,3)][1,]
#   means[i,3:6] <- colMeans(analysis.out[seq(from=i, to=50, by=5),4:7])
# }
# for(i in 1:5){
#   means[i+5,1:2] <- analysis.out[seq(from=i+50, to=100, by=5),c(1,3)][1,]
#   means[i+5,3:6] <- colMeans(analysis.out[seq(from=i+50, to=100, by=5),4:7])
# }
# for(i in 1:5){
#   means[i+10,1:2] <- analysis.out[seq(from=i+100, to=150, by=5),c(1,3)][1,]
#   means[i+10,3:6] <- colMeans(analysis.out[seq(from=i+100, to=150, by=5),4:7])
# }
# for(i in 1:5){
#   means[i+15,1:2] <- analysis.out[seq(from=i+150, to=200, by=5),c(1,3)][1,]
#   means[i+15,3:6] <- colMeans(analysis.out[seq(from=i+150, to=200, by=5),4:7])
# }
# for(i in 1:5){
#   means[i+20,1:2] <- analysis.out[seq(from=i+201, to=250, by=5),c(1,3)][1,]
#   means[i+20,3:6] <- colMeans(analysis.out[seq(from=i+201, to=250, by=5),4:7])
# }
# for(i in 1:5){
#   means[i+25,1:2] <- analysis.out[seq(from=i+251, to=300, by=5),c(1,3)][1,]
#   means[i+25,3:6] <- colMeans(analysis.out[seq(from=i+251, to=300, by=5),4:7])
# }
# 
# 
# ald.row <- which(means$method == "ald")
# ald2.row <- which(means$method == "ald2")
# ald5.row <- which(means$method == "ald5")
# des.row <- which(means$method == "des")
# des5.row <- which(means$method == "des5")
