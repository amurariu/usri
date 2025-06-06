lim.fun <- function(y, data, conditions, nloop=100){
  conditions_p <- conditions
  conds <- as.data.frame(conditions_p)
  
  thin.data.out <- list() 
  data.out.limma.u <- list() 
  data.out.limma.r <- list() 
  data.out.limma.p <- list() 
  
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
    
    #randomized without FP addition PD1
    #y <- calcNormFactors(y) #pd_1/immuno.data
    logCPM <- cpm (y_pd1, log = TRUE)
    design<-model.matrix(~condsp)
    v<-voom(y,design)
    fit <- lmFit(v,design)
    fit <- eBayes(fit)
    res.lim <- topTable(fit)
    
    #randomized and FP addition PD1
    data.p <- as.data.frame(thin$mat)
    yfd <- calcNormFactors(data.p) #pd_1/immuno.data
    design<-model.matrix(~condsp)
    v<-voom(data.p,design)
    fit <- lmFit(v,design)
    fit <- eBayes(fit)
    res.lim <- topTable(fit)
    
  }
  print("done loop")
  
  #unpermuted PD1
  #randomized without FP addition PD1
  y <- calcNormFactors(y) #pd_1/immuno.data
  design<-model.matrix(~conds) #normal conditions
  v<-voom(y,design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  res.lim <- topTable(fit)
 
  
  return(list(conditions=conditions_p, thin.data=thin.data.out, data.out.limma.r = res.lim))
}
