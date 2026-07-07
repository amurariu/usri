# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald2.fun <- function(data, conditions, nloop=100, gamma, prop_null, mean, nsample=128, std=2){

  thin.data.out.aldex <- list() 
  data.out.aldex.r <- list()
  data.out.aldex.t <- list() 
  
  #for loop
  for (i in 1:nloop){
    seed = 20 + i
    set.seed(seed)
    
    message(paste0("\nLoop iteration: ", i))
    
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    message('Thinning data...')
    thin <- thin_2group(data, prop_null=prop_null, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = mean, sd = std))
    thin.data.out.aldex[[i]] <- thin
    conds_th <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    data_th <- thin$mat
    
    #randomized without FP addition
    message("Running ALDEx2 (scale = ", gamma, ") on original data with randomised groups...")
    aldex.r <- aldex(data, conditions=conds_th, gamma = gamma) #uses original dataset but permuted conditions
    data.out.aldex.r[[i]] <- aldex.r
    
    #randomized with FP addition 
    message("Running ALDEx2 (scale = ", gamma, ") on thinned data with randomised groups...")
    aldex.t <- aldex(data_th, conditions=conds_th, gamma = gamma, mc.samples = nsample) #uses new dataset with permuted conditions
    data.out.aldex.t[[i]] <- aldex.t
  }
  
  message("\nAll ", i, " iterations of loop finished\n")
  message("Running ALDEx2 (scale = ", gamma, ") on original data with non-randomised groups...")
  
  set.seed(2025)
  #unpermuted
  aldex.u <- aldex(data, conditions=conditions, gamma = gamma)

  return(list(conditions=conditions, r.data=data.out.aldex.r, thin.data=thin.data.out.aldex, u.data=aldex.u, t.data=data.out.aldex.t))
}

