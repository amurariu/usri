lim.fun <- function(data, conditions, nloop=100){
  
  thin.data.out <- list() 
  data.out.limma.r <- list() 
  data.out.limma.t <- list() 
  
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
    conds_th <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    data_th <- thin$mat
    
    #randomized without FP addition
    designr<-model.matrix(~conds_th)
    vr<-voom(data,designr) 
    fitr <- lmFit(vr,designr)
    fitr <- eBayes(fitr)
    data.out.limma.r[[i]] <- topTable(fitr)
    
    #randomized and FP addition 
    designt<-model.matrix(~conds_th)
    vt<-voom(data_th,designt) #thinned data+conditions
    fitt <- lmFit(vt,designt)
    fitt <- eBayes(fitt)
    data.out.limma.t[[i]] <- topTable(fitt, coef = ncol(designt))
    
  }
  print("done loop")
  
  #unpermuted PD1 
  set.seed(2025)
  designu<-model.matrix(~conditions)
  vu<-voom(data,designu) #original conditions + original data
  fitu <- lmFit(vu,designu)
  fitu <- eBayes(fitu)
  res.lim.u <- topTable(fitu, coef = ncol(designu))
  
  return(list(conditions=conditions, thin.data=thin.data.out, data.out.limma.r = data.out.limma.r, data.out.limma.t = data.out.limma.t, data.out.limma.u = res.lim.u))
}
