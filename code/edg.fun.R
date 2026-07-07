# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
edg.fun <- function(data, conditions, nloop = 100, prop_null = 0.95, mean = 0, std = 2){
  
  thin.data.out.edger <- list() #change name of list here-----------
  data.out.edgeR.u <- list() 
  data.out.edgeR.r <- list() 
  data.out.edgeR.t <- list() 
  
  #for loop
  for (i in 1:nloop){
    
    message(paste0("\nLoop iteration: ", i))
    
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    # ensure that thinning is the different for all iterations
    # use in each function, so that it is consistent dataset for each tool
    seed = 20 + i
    set.seed(seed)
    
    message('Thinning data...')
    thin <- thin_2group(data, prop_null = prop_null, alpha=0,
                               signal_fun = stats::rnorm, 
                               signal_params = list(mean = mean, sd = std))
    thin.data.out.edger[[i]] <- thin
    conds_th <- as.vector(thin$designmat)   #thinned and randomized conditions
    data_th <- thin$mat #thinned and randomized data
    
    #edgeR analysis
    group_th <- factor(conds_th)
    design_th <- model.matrix(~group_th) #use group randomization from seqgendiff
    
    #randomized without FP addition
    message("Running edgeR on original data with randomised groups...")
    fit_r <- glmQLFit(data,design_th) #uses original data (ie. no TP added)
    qlf_r <- glmQLFTest(fit_r,coef=2)
    edg.r<-topTags(qlf_r, n=nrow(data), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    data.out.edgeR.r[[i]] <- as.data.frame(edg.r[[1]])
    
    #randomized with FP addition 
    message("Running edgeR on thinned data with randomised groups...")
    fit_t <- glmQLFit(data_th,design_th)
    qlf_t <- glmQLFTest(fit_t,coef=2)
    edg.t<-topTags(qlf_t, n=nrow(data_th), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    data.out.edgeR.t[[i]] <- as.data.frame(edg.t[[1]])
  }
  
  message("\nAll ", i, " iterations of loop finished\n")
  
  set.seed(2025)
  
  #unpermuted 
  message("Running edgeR on original data with non-randomised groups...")
  group_u<-factor(conditions)
  design_u <- model.matrix(~group_u)
  fit_u <- glmQLFit(data,design_u)
  qlf_u <- glmQLFTest(fit_u,coef=2)
  edg.u<-topTags(qlf_u, n=nrow(data), adjust.method = "BH", sort.by = "none", p.value = 1) 
  
  data.out.edgeR.u <- as.data.frame(edg.u[[1]])
  
  return(list(conditions=conditions, thin.data=thin.data.out.edger, u.data=data.out.edgeR.u, r.data=data.out.edgeR.r, t.data=data.out.edgeR.t))
}
