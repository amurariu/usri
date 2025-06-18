lim.fun <- function(data, conditions, nloop=100){
  
  condition <- conditions
  conds <- as.vector(condition)
  
  thin.data.out <- list() 
  data.out.limma.r <- list() 
  data.out.limma.p <- list() 
  
  #for loop
  for (i in 1:nloop){
    print(i)
    seed = 20 + i
    set.seed(seed)
    
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    thin <- thin_2group(data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out[[i]] <- thin
    condsp <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    datasp <- thin$mat
    
    #randomized without FP addition
    designp<-model.matrix(~condsp)
    vp<-voom(data,designp) 
    fitp <- lmFit(vp,designp)
    fitp <- eBayes(fitp)
    data.out.limma.p[[i]] <- topTable(fitp)
    
    #randomized and FP addition 
    designr<-model.matrix(~condsp)
    vr<-voom(datasp,designr) #thinned data+conditions
    fitr <- lmFit(vr,designr)
    fitr <- eBayes(fitr)
    data.out.limma.r[[i]] <- topTable(fitr, coef = ncol(designr))
    
  }
  print("done loop")
  
  #unpermuted PD1 
  set.seed(2025)
  designu<-model.matrix(~conds)
  vu<-voom(data,designu) #original conditions + original data
  fitu <- lmFit(vu,designu)
  fitu <- eBayes(fitu)
  res.lim.u <- topTable(fitu, coef = ncol(designu))
  
  return(list(conditions=conditions, thin.data=thin.data.out, data.out.limma.r = data.out.limma.r, data.out.limma.p = data.out.limma.p, data.out.limma.u = res.lim.u))
}
