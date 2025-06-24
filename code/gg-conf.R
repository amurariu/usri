
# tests
# fp.DES <- get_confusion(immuno.data.DESeq, prog="DESeq")
# fp.edge <- get_confusion(immuno.data.edgeR, prog="edgeR")
# fp.limma <- get_confusion(immuno.data.limma, prog="limma")
# fp.A0 <- get_confusion(immuno.data_0.aldex2, prog="ALDEx2") # t-test
# fp.A0.w <- get_confusion(immuno.data_0.aldex2, prog="ALDEx2.wi") # wilcox
# can check the output with summary(output)

get_confusion <- function(input=NULL, prog=NULL){
	#if(input==NULL){stop("input dataset")}
	#if(type==NULL){stop("input dataset type")}
	if(prog=="DESeq" || prog=="deseq"){
		lfc = 2 # fold change column
		padj = 6 # BH adjusted p column
	}else if(prog=="edgeR" || prog == "edger"){
		lfc = 1 # fold change column
		padj = 5 # BH adjusted p column
	}else if(prog=="limma"){
		lfc = 1
		padj = 5
	}else if(prog=="ALDEx2" || prog=="aldex" || prog=="aldex2"){
		lfc=4
		padj=9
	}else if(prog=="ALDEx2.wi"){
		lfc=4
		padj=11
	}
	
	# get the confusion matrix for the p.data

	coeff <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
	# count of all FP in random
	TPFPR.mat = matrix(data=NA, nrow=length(coeff), ncol=4)
	colnames(TPFPR.mat) <- c("TPR0", "FDR0", "TPR5","FDR5")
	rownames(TPFPR.mat) <- coeff
	
	data.out <- list()
	
	for(j in 1:length(coeff)){ # coeffient counter
		# for each coefficient
		raw.counts <- matrix(data=NA, nrow=length(input$t.data), ncol=7)
		diff.counts <- matrix(data=NA, nrow=length(input$t.data), ncol=7) # same with difference of 0.5 fold
		colnames(raw.counts) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(diff.counts) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
for(i in 1:length(input$t.data)){ # replicate counter
		# all positives in the random data are FP		
		TP.coeff <- which(abs(input$thin.data[[i]]$coefmat) >= coeff[j]) # true null
		TN.coeff <- which(abs(input$thin.data[[i]]$coefmat) == 0 )
		# model
		raw.counts[i,5] <- length(TP.coeff)
		diff.counts[i,5] <- length(TP.coeff)
		# TP
		raw.counts[i,1] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05), TP.coeff))
		diff.counts[i,1] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05 & 
		  abs(input$t.data[[i]][,lfc]) > 0.5), TP.coeff))
		# FP
		raw.counts[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05), TN.coeff))
		diff.counts[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05 & 
		  abs(input$t.data[[i]][,lfc]) > 0.5), TN.coeff))
		# TN
		raw.counts[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05), TN.coeff))
		diff.counts[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05 & 
		  abs(input$t.data[[i]][,lfc]) < 0.5), TN.coeff))
		# FN
		raw.counts[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05), TP.coeff))
		diff.counts[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05 & 
		  abs(input$t.data[[i]][,lfc]) < 0.5), TP.coeff))
		# TPR
		raw.counts[i,6] <- round(raw.counts[i,1]/raw.counts[i,5],3)
		diff.counts[i,6] <- round(diff.counts[i,1]/diff.counts[i,5],3)
		# FDR
		raw.counts[i,7] <- round(raw.counts[i,2]/(raw.counts[i,1] + raw.counts[i,2]),3)
		diff.counts[i,7] <- round(diff.counts[i,2]/(diff.counts[i,1] + diff.counts[i,2]),3)
 		}
 		data.out[[j]] <- raw.counts
 		
 		TPFPR.mat[j,1] <- mean(raw.counts[,6])
 		TPFPR.mat[j,2] <- mean(raw.counts[,7])
 		TPFPR.mat[j,3] <- mean(diff.counts[,6])
 		TPFPR.mat[j,4] <- mean(diff.counts[,7])
 		
	}
	data.out[[length(coeff)+1]] <- TPFPR.mat
	names(data.out) <- c(coeff, "TPFPR")
	return(data.out)
	
}