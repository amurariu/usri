# example function for DESeq2
# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
des.fun <- function(data, conditions, nloop=4){
    
	#assign(paste("perf.a", "1", sep=""),5)
    #perf.a1
	conditions_p <- conditions
	conds <- data.frame(conditions_p)
	
	thin.data.out <- list() 
	data.out.deseq.u <- list() 
	data.out.deseq.r <- list() 
	data.out.deseq.p <- list() 
	
	#for loop
	for (i in 1:nloop){
	  print(i)
	  #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
	  #generate thin_2group for each dataset as well as labelling for conditions and new dataset
	  
	  #PD1
	  thin <- thin_2group(data, prop_null=0.95, alpha=0,
	                      signal_fun = stats::rnorm, 
	                      signal_params = list(mean = 0, sd = 2))
	  thin.data.out[[i]] <- thin
	  condsp <- as.vector(thin$designmat)   # permuted and thinned conditions and data
	  datasp <- thin$mat
	  
	  #DESeq2 analysis
	  #randomized without FP addition PD1
	  dds.rp.deseq  <- DESeqDataSetFromMatrix(countData = data,  #uses original data (no TP added)
									   colData = data.frame(condsp), #uses data randomization order from thin
									   design = ~ condsp)
	  dds.rp.deseq <- DESeq(dds.rp.deseq, quiet=T)
	  res.rp.deseq <- results(dds.rp.deseq)
	  data.out.deseq.r[[i]] <- as.data.frame(res.rp.deseq@listData) #added [[i]] and referenced list in line prior
	  
	  #randomized with FP addition PD1
	  dds.thp.deseq  <- DESeqDataSetFromMatrix(countData = datasp,
										colData = data.frame(condsp),
										design = ~ condsp)
	  dds.thp.deseq <- DESeq(dds.thp.deseq, quiet=T)
	  res.thp.deseq <- results(dds.thp.deseq)
	  data.out.deseq.p[[i]] <- as.data.frame(res.thp.deseq@listData)
	}
	print("done loop")

	
	#unpermuted PD1
	dds.up.deseq  <- DESeqDataSetFromMatrix(countData = data,
									 colData = conds,
									 design = ~ conditions_p)
	dds.up.deseq <- DESeq(dds.up.deseq, quiet=T)
	data.out.deseq.u <- results(dds.up.deseq)
	data.des.u <- as.data.frame(data.out.deseq.u@listData)
	
	return(list(conditions=conditions_p, thin.data=thin.data.out, u.data=data.des.u, r.data=data.out.deseq.r, p.data=data.out.deseq.p))
}
