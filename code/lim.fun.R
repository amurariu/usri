lim.fun <- function(y, data, conditions, nloop=100){
  conditions_p <- conditions
  conds <- as.data.frame(conditions_p, row.names = NULL, optional = FALSE,
                         make.names = TRUE,
                         stringsAsFactors = FALSE)
  
  thin.data.out <- list() 
  data.out.deseq.u <- list() 
  data.out.deseq.r <- list() 
  data.out.deseq.p <- list() 
  
  #for loop
  for (i in 1:nloop){
    print(i)
    seed = 20 + i
    set.seed(seed)
    
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin <- thin_2group(data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out[[i]] <- thin
    condsp <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    datasp <- thin$mat
    
    #limma analysis
    #randomized without FP addition PD1
  
    y <- calcNormFactors(y)
    design<-model.matrix(~condsp)
    v<-voom(y,design)
    fit <- lmFit(v,design)
    fit <- eBayes(fit)
    res.lim = topTable(fit)
    
  }
  print("done loop")
  
  
  #unpermuted PD1
 
  
  return(list(conditions=conditions_p, thin.data=thin.data.out, u.data=data.des.u, r.data=data.out.deseq.r, p.data=data.out.deseq.p))
}
