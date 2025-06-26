### 
# outputs a list of the actual data for the TPR/FDR calculations
# and a summary of that data 
# two definitions of negatives
# first, is negative if modelled LFC (coeff) is less than a prescribed LFC
# second, is negative if there is no modelled LFC
###
# example calling
# fp.DES <- get_confusion(immuno.data.DESeq, prog="DESeq")
# fp.edge <- get_confusion(immuno.data.edgeR, prog="edgeR")
# fp.limma <- get_confusion(immuno.data.limma, prog="limma")
# fp.A0 <- get_confusion(immuno.data_0.aldex2, prog="ALDEx2") # t-test
# fp.A0.w <- get_confusion(immuno.data_0.aldex2, prog="ALDEx2.wi") # wilcox
# can check the output with summary(output$TPFPR)

###

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
	
	# get the confusion matrix for the t.data

	coeff <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
	# count of all FP in random
	TPFPR.mat = matrix(data=NA, nrow=length(coeff), ncol=8)
	colnames(TPFPR.mat) <- c("cTPR0", "cFDR0", "cTPR5","cFDR5","zTPR0", "zFDR0", "zTPR5","zFDR5")
	rownames(TPFPR.mat) <- coeff
	
	raw_coeff <- list()
	diff_coeff <- list()
	raw_zero <- list()
	diff_zero <- list()
	
	for(j in 1:length(coeff)){ # coeffient counter
		# for each coefficient
		raw.coeff <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7))
		diff.coeff <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7)) # same with difference of 0.5 fold
		raw.zero <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7))
		diff.zero <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7)) # same with difference of 0.5 fold
		colnames(raw.coeff) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(diff.coeff) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(raw.zero) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(diff.zero) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
			
		for(i in 1:length(input$t.data)){ # replicate counter
			# all positives in the random data are FP		
			TP.coeff <- which(abs(input$thin.data[[i]]$coefmat) >= coeff[j]) # true null
			TN.coeff <- which(abs(input$thin.data[[i]]$coefmat) < coeff[j])
			TN.zero <- which(abs(input$thin.data[[i]]$coefmat) == 0)
			# model is the same for all 
			raw.coeff[i,5] <- length(TP.coeff)
			diff.coeff[i,5] <- length(TP.coeff)
			raw.zero[i,5] <- length(TP.coeff)
			diff.zero[i,5] <- length(TP.coeff)
			
			# TP, the same for all
			raw.coeff[i,1] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05), TP.coeff))
			raw.zero[i,1] <- raw.coeff[i,1]
			diff.coeff[i,1] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05 & 
			  abs(input$t.data[[i]][,lfc]) > 0.5), TP.coeff))
			diff.zero[i,1] <- diff.coeff[i,1]
			
			# FP, differ by negative definition
			raw.coeff[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05), TN.coeff))
			raw.zero[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05), TN.zero))
			diff.coeff[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05 & 
			                                            abs(input$t.data[[i]][,lfc]) > 0.5), TN.coeff))
			diff.zero[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < 0.05 & 
			  abs(input$t.data[[i]][,lfc]) > 0.5), TN.zero))
			
			# TN, differ by negative definition
			raw.coeff[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05), TN.coeff))
			diff.coeff[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05 & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TN.coeff))
			raw.zero[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05), TN.zero))
			diff.zero[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05 & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TN.zero))
			
			# FN, same for all
			raw.coeff[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05), TP.coeff))
			diff.coeff[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05 & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TP.coeff))
			raw.zero[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05), TP.coeff))
			diff.zero[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > 0.05 & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TP.coeff))
			# TPR
			raw.coeff[i,6] <- round(raw.coeff[i,1]/raw.coeff[i,5],3)
			diff.coeff[i,6] <- round(diff.coeff[i,1]/diff.coeff[i,5],3)
			raw.zero[i,6] <- round(raw.zero[i,1]/raw.zero[i,5],3)
			diff.zero[i,6] <- round(diff.zero[i,1]/diff.zero[i,5],3)
			# FDR
			raw.coeff[i,7] <- round(raw.coeff[i,2]/(raw.coeff[i,1] + raw.coeff[i,2]),3)
			diff.coeff[i,7] <- round(diff.coeff[i,2]/(diff.coeff[i,1] + diff.coeff[i,2]),3)
			raw.zero[i,7] <- round(raw.zero[i,2]/(raw.zero[i,1] + raw.zero[i,2]),3)
			diff.zero[i,7] <- round(diff.zero[i,2]/(diff.zero[i,1] + diff.zero[i,2]),3)
			}
			raw_coeff[[j]] <- raw.coeff
			diff_coeff[[j]] <- diff.coeff
			raw_zero[[j]] <- raw.zero
			diff_zero[[j]] <- diff.coeff
			
			TPFPR.mat[j,1] <- mean(raw.coeff[,6]) # TPR
			TPFPR.mat[j,2] <- mean(raw.coeff[,7]) # FDR
			TPFPR.mat[j,3] <- mean(diff.coeff[,6])
			TPFPR.mat[j,4] <- mean(diff.coeff[,7])
			TPFPR.mat[j,5] <- mean(raw.zero[,6]) # TPR
			TPFPR.mat[j,6] <- mean(raw.zero[,7]) # FDR
			TPFPR.mat[j,7] <- mean(diff.zero[,6])
			TPFPR.mat[j,8] <- mean(diff.zero[,7])
			
		}
		names(raw_coeff) <- coeff
		names(diff_coeff) <- coeff
		names(raw_zero) <- coeff
		names(diff_zero) <- coeff
		data.out <- list(raw_coeff, diff_coeff, raw_zero, diff_zero, TPFPR.mat)
		names(data.out) <- c("raw_coeff","diff_coeff","raw_zero","diff_zero","TPFPR")
		return(data.out)
}