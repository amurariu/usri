# function for generating DESeq2 data x100

# Author: Andreea Murariu
# Last edited: 2025-06-16

# data: matrix of raw counts to be analysed, columns = samples, rows = feature
# conditions: original grouping variable vector
# nloop: number of randomisation instances for randomised data

des.fun <- function(data, nloop=100, conditions){
	
	thin.data.out <- list() 
	data.out.deseq.u <- list()
	data.out.deseq.t <- list() 
	data.out.deseq.r <- list() 
	
	# initialise loop
	for (i in 1:nloop){
	  
	  message(paste0("Loop iteration: ", i))
	  
	  # thin_2group adds rnorm noise to 5% of the transcripts, generating 
	  # 'true positives' in the data, using an approach called binomial thinning.
	  # Also permutes the grouping labels so that samples are randomly placed in
	  # one of two groups.
	  seed = 20 + i
	  set.seed(seed)
	  
	  # binomial thinning of data
	  thin <- thin_2group(data, prop_null=0.95, alpha=0,
	                      signal_fun = stats::rnorm, 
	                      signal_params = list(mean = 0, sd = 2))
	  
	  thin.data.out[[i]] <- thin
	  
	  # extract randomised conditions matrix and convert to factor for DESeq2
	  # (factor required or else deseq will specify a model with increasing fold
	  # change for higher values)
	  condsr <- as.data.frame(thin$designmat)
	  condsr[,1] <- as.factor(condsr[,1])
	  
	  # extract thinned read count data (i.e. contains 'true' positives)
	  datasp <- thin$mat
	  
	  # DESeq2 analysis
	  #------- original data with randomised groups
	  message("Running DESeq2 on original data with randomised groups...")
	  
	  # pull colname of randomised conditions matrix as formula
	  form.r <- as.formula(paste0("~", colnames(condsr)))
	  
	  # run DESeq2
	  dds.r.deseq  <- DESeqDataSetFromMatrix(countData = data, design = form.r,
	                                          colData = data.frame(condsr))
	  
	  dds.r.deseq <- DESeq(dds.r.deseq, quiet=T)
	  res.r.deseq <- results(dds.r.deseq)
	  data.out.deseq.r[[i]] <- as.data.frame(res.r.deseq@listData)
	  
	  #------- thinned data with randomised groups
	  message("Running DESeq2 on thinned data with randomised groups...")
	  
	  dds.th.deseq  <- DESeqDataSetFromMatrix(countData = datasp, design = form.r,
	                                           colData = data.frame(condsr))
	  
	  dds.th.deseq <- DESeq(dds.th.deseq, quiet=T)
	  res.th.deseq <- results(dds.th.deseq)
	  data.out.deseq.t[[i]] <- as.data.frame(res.th.deseq@listData)
	}
	
	message("\nAll ", i, " iterations of loop finished\n")
	message("Running DESeq2 on original data with non-randomised groups...")
	
	# ensure conditions is a data frame and convert data to factor
	conds <- as.data.frame(conditions)
	conds[,1] <- as.factor(conds[,1])
	
	#------- original data with non-randomised groupings
	
	# pull colname of original conditions matrix as formula
	form <- as.formula(paste0("~", colnames(conds)))
	
	# run DESeq2
	set.seed(2025)
	dds.u.deseq  <- DESeqDataSetFromMatrix(countData = data, colData = conds,
	                                        design = form)
	
	dds.u.deseq <- DESeq(dds.u.deseq, quiet=T)
	data.des.u <- results(dds.u.deseq)
	data.out.deseq.u <- as.data.frame(data.des.u@listData)
	
	# return a list containing the following:
	#   - original conditions vector
	#   - list of 100x thinned data (itself a ThinData object)
	#   - list of 1x DESeq2 results: original data with non-randomised groups
	#   - list of 100x DESeq2 results: original data with randomised groups
	#   - list of 100x DESeq2 results: thinned data with randomised groups
	return(list(conditions=conditions, thin.data=thin.data.out,
	            u.data=data.out.deseq.u, r.data=data.out.deseq.r, t.data=data.out.deseq.t))
}
