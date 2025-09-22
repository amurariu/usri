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

get_confusion <- function(input=NULL, prog=NULL, FDR=0.1){
	#if(input==NULL){stop("input dataset")}
	#if(type==NULL){stop("input dataset type")}
	
  # pull indices of columns containing log fold-change (or similar measure of 
  # difference between groups) and adjusted P-value
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
	}else if (prog=="ALDEx3" || prog=="aldex3") {
	  lfc= 3
	  padj=5
	}
	
	# get the confusion matrix for the t.data

  # make vector of prescribed log-fold change thresholds (i.e. which represent 
  # modelled log fold-change injected into data from binomial thinning using the
  # 'seqgendiff' package)  
	coeff <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) #modelled LFC from thinning
	
	# make a matrix to hold counts of all the false-positives and true-positives
	# count of all FP in random
	TPFPR.mat = matrix(data=NA, nrow=length(coeff), ncol=8)
	colnames(TPFPR.mat) <- c("cTPR0", "cFDR0", "cTPR5","cFDR5","zTPR0", "zFDR0", "zTPR5","zFDR5")
	rownames(TPFPR.mat) <- coeff
	
	# make 4 lists to hold true/false positive/negative features of 100 loop 
	# iterations for thinned data
	raw_coeff <- list()   # explanation
	diff_coeff <- list()  # explanation
	raw_zero <- list()    # explanation
	diff_zero <- list()   # explanation
	
	for(j in 1:length(coeff)){
		
	  # make dataframes & set column names to hold following data for the current
	  # coefficient threshold:
	  raw.coeff <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7))
		diff.coeff <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7)) # same with difference of 0.5 fold
		raw.zero <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7))
		diff.zero <- as.data.frame(matrix(data=NA, nrow=length(input$t.data), ncol=7)) # same with difference of 0.5 fold
		
		colnames(raw.coeff) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(diff.coeff) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(raw.zero) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
		colnames(diff.zero) <- c("TP", "FP", "TN", "FN","model","TPR","FDR")
			
		for(i in 1:length(input$t.data)){ # replicate counter
			
		  # all positives in the thinned data are false-positives as the sample
		  # groupings have been randomised; however we can distinguish TPs and TNs
		  # for each coefficient threshold. Two kinds of TN features: those not
		  # altered by thinning (zeros in coefmat) and those whose coefficient
		  # was altered by an amount less than the coefficient threshold (for the
		  # current loop)
			TP.coeff <- which(abs(input$thin.data[[i]]$coefmat) >= coeff[j]) # true positives; FC altered by thinning & >threshold
			TN.coeff <- which(abs(input$thin.data[[i]]$coefmat) < coeff[j]) # true negatives; FC altered by thinning & <threshold AND FC not altered
			TN.zero <- which(abs(input$thin.data[[i]]$coefmat) == 0) # true negatives; FC not altered by thinning
		
			# 'model' is the same for all: number of features where the MODELLED
			# coefficient is > coefficient threshold for current loop
			raw.coeff[i,5] <- length(TP.coeff)
			diff.coeff[i,5] <- length(TP.coeff)
			raw.zero[i,5] <- length(TP.coeff)
			diff.zero[i,5] <- length(TP.coeff)
			
			# true positives in the MEASURED data (thinned & randomised) are identical 
			# for both coeff/zero but will differ between raw/diff based on an 
			# arbitrary LFC threshold of >0.5
			
			# raw.coeff/raw.zero:   features with a MODELLED FC altered by thinning
			#                       that is >threshold, with a MEASURED P value <FDR
			
			# diff.coeff/diff.zero: features with a MODELLED FC altered by thinning
			#                       that is >threshold, with a MEASURED P value <FDR
			#                       and a MEASURED log fold-change >0.5
			
			raw.coeff[i,1] <- length(intersect(which(input$t.data[[i]][,padj] < FDR), TP.coeff))
			raw.zero[i,1] <- raw.coeff[i,1]
			diff.coeff[i,1] <- length(intersect(which(input$t.data[[i]][,padj] < FDR & 
			  abs(input$t.data[[i]][,lfc]) > 0.5), TP.coeff))
			diff.zero[i,1] <- diff.coeff[i,1]
			
			# false positives in the MEASURED data (thinned & randomised) will differ 
			# between coeff/zero based on how negatives are defined, and between
			# raw/diff based on an arbitrary LFC threshold of >0.5
			
			# raw.coeff: features with a MODELLED FC altered by thinning that is
			#            <threshold (as well as FC not altered), with a MEASURED P value <FDR
			
			# raw.zero:  features whose fold changes were NOT altered by thinning,
			#            but whose MEASURED P values are <FDR
			
			# diff.coeff/diff.zero: same as respective definitions above, but all FPs
			#                       must also have a MEASURED log fold-change >0.5
			
			raw.coeff[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < FDR), TN.coeff))
			raw.zero[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < FDR), TN.zero))
			diff.coeff[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < FDR & 
			                                            abs(input$t.data[[i]][,lfc]) > 0.5), TN.coeff))
			diff.zero[i,2] <- length(intersect(which(input$t.data[[i]][,padj] < FDR & 
			                                           abs(input$t.data[[i]][,lfc]) > 0.5), TN.zero))
			
			# true negatives in the MEASURED data (thinned & randomised) will differ
			# between coeff/zero based on how negatives are defined, and between
			# raw/diff based on an arbitrary LFC threshold of >0.5
			
			# raw.coeff: features with a MODELLED FC altered by thinning that is
			#            <threshold (as well as FC not altered), with a MEASURED P value >FDR
			
			# raw.zero:  features with a MODELLED FC that was NOT altered by thinning 
			#            and with a MEASURED P value >FDR
			
			# diff.coeff/diff.zero: same as respective definitions above, but all TNs
			#                       must also have a MEASURED log fold-change <0.5
			
			raw.coeff[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > FDR), TN.coeff))
			diff.coeff[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > FDR & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TN.coeff))
			raw.zero[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > FDR), TN.zero))
			diff.zero[i,3] <- length(intersect(which(input$t.data[[i]][,padj] > FDR & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TN.zero))
			
			# false negatives in the MEASURED data (thinned & randomised) will be
			# identical for both coeff/zero but will differ between raw/diff based on 
			# an arbitrary LFC threshold of <0.5
			
			# raw.coeff/raw.zero:   features with a MODELLED FC altered by thinning
			#                       that is >threshold, with a MEASURED P value >FDR
			
			# diff.coeff/diff.zero: features with a MODELLED FC altered by thinning
			#                       that is >threshold, with a MEASURED P value >FDR
			#                       and a MEASURED log fold-change >0.5
			
			raw.coeff[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > FDR), TP.coeff))
			diff.coeff[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > FDR & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TP.coeff))
			raw.zero[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > FDR), TP.coeff))
			diff.zero[i,4] <- length(intersect(which(input$t.data[[i]][,padj] > FDR & 
			  abs(input$t.data[[i]][,lfc]) < 0.5), TP.coeff))
		
			# the true positive rate (TPR) is the number of of true positives in the
			# MEASURED data divided by the number of features in the MODELLED data
			# which have a coefficient greater than the prescribed threshold (i.e how
			# many known TPs are we detecting)
			
			raw.coeff[i,6] <- round(raw.coeff[i,1]/raw.coeff[i,5],3)
			diff.coeff[i,6] <- round(diff.coeff[i,1]/diff.coeff[i,5],3)
			raw.zero[i,6] <- round(raw.zero[i,1]/raw.zero[i,5],3)
			diff.zero[i,6] <- round(diff.zero[i,1]/diff.zero[i,5],3)
			
			# the false discovery rate (FDR) is the number of false positives in the
			# MEASURED data divided by number of true positives and false positives
			# in the MEASURED data (i.e. how many times are we calling something truly
			# different between groups when we know it is not)
			
			raw.coeff[i,7] <- round(raw.coeff[i,2]/(raw.coeff[i,1] + raw.coeff[i,2]),3)
			diff.coeff[i,7] <- round(diff.coeff[i,2]/(diff.coeff[i,1] + diff.coeff[i,2]),3)
			raw.zero[i,7] <- round(raw.zero[i,2]/(raw.zero[i,1] + raw.zero[i,2]),3)
			diff.zero[i,7] <- round(diff.zero[i,2]/(diff.zero[i,1] + diff.zero[i,2]),3)
			}
		
		# add all 100 iterations of the TP/FP/TN/FN/modelCoeffs/TPR/FDR to their
		# respective lists
		raw_coeff[[j]] <- raw.coeff
		diff_coeff[[j]] <- diff.coeff
		raw_zero[[j]] <- raw.zero
		diff_zero[[j]] <- diff.coeff
			
		# calculate the average of the true positive rate and false discovery rate
		# across all 100 instances and store in matrix
		TPFPR.mat[j,1] <- mean(raw.coeff[,6])   # cTPR0
		TPFPR.mat[j,2] <- mean(raw.coeff[,7])   # cFDR0
		
		TPFPR.mat[j,3] <- mean(diff.coeff[,6])  # cTPR5
		TPFPR.mat[j,4] <- mean(diff.coeff[,7])  # cFDR5
		
		TPFPR.mat[j,5] <- mean(raw.zero[,6])    # zTPR0
		TPFPR.mat[j,6] <- mean(raw.zero[,7])    # zFDR0
		
		TPFPR.mat[j,7] <- mean(diff.zero[,6])   # zTPR5
		TPFPR.mat[j,8] <- mean(diff.zero[,7])   # zFDR5
		}
	
	# set names of the list elements to the prescribed coefficient thresholds
	names(raw_coeff) <- coeff
	names(diff_coeff) <- coeff
	names(raw_zero) <- coeff
	names(diff_zero) <- coeff
	
	# build final list of data, set names and return object
	data.out <- list(raw_coeff, diff_coeff, raw_zero, diff_zero, TPFPR.mat)
	names(data.out) <- c("raw_coeff","diff_coeff","raw_zero","diff_zero","TPFPR")
	return(data.out)
}