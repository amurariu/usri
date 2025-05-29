# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald2.fun <- function(data, conditions, nloop=100, gamma){
  set.seed(4)
  #assign(paste("perf.a", "1", sep=""),5)
  #perf.a1
  #gam <- c(1e-3, 0.2, 0.5) #contains different scale values

  thin.data.out.aldex <- list() #change name of list here-----------
  data.out.aldex.u <- list() 
  data.out.aldex.r <- list() 
  data.out.aldex.p <- list() 
  
  #for loop
  for (i in 1:nloop){
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin <- thin_2group(data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out.aldex[[i]] <- thin
    condsp <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    datasp <- thin$mat
    
    #randomized without FP addition PD1
    xrp.aldex <- aldex(data, conditions=condsp, gamma = gamma) #uses original dataset but permuted conditions
    data.out.aldex.r[[i]] <- xrp.aldex
    
    #randomized with FP addition PD1
    xpp.aldex <- aldex(datasp, conditions=condsp, gamma = gamma) #uses new dataset with permuted conditions
    data.out.aldex.p[[i]] <- xpp.aldex
  }
  
  print("done loop")
  
  #unpermuted PD1
  xup.aldex <- aldex(data, conditions=conditions, gamma = gamma)

  return(list(conditions=conditions, thin.data=thin.data.out.aldex, u.data=xup.aldex, r.data=data.out.aldex.r, p.data=data.out.aldex.p))
}
