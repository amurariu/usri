# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald2.fun <- function(data, conditions, nloop=100, gamma){

  thin.data.out.aldex <- list() 
  data.out.aldex.r <- list() 
  data.out.aldex.t <- list() 
  
  #for loop
  for (i in 1:nloop){
    seed = 20 + i
    set.seed(seed)
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin <- thin_2group(data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out.aldex[[i]] <- thin
    conds_th <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    data_th <- thin$mat
    
    #randomized without FP addition PD1
    aldex.r <- aldex(data, conditions=conds_th, gamma = gamma) #uses original dataset but permuted conditions
    data.out.aldex.r[[i]] <- aldex.r
    
    #randomized with FP addition PD1
    aldex.t <- aldex(data_th, conditions=conds_th, gamma = gamma) #uses new dataset with permuted conditions
    data.out.aldex.t[[i]] <- aldex.t
  }
  
  print("done loop")
  
  #unpermuted PD1
  aldex.u <- aldex(data, conditions=conditions, gamma = gamma)

  return(list(conditions=conditions, thin.data=thin.data.out.aldex, u.data=aldex.u, r.data=data.out.aldex.r, t.data=data.out.aldex.t))
}
