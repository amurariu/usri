# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald3.fun <- function(data, conditions, nloop=100, gamma){
  set.seed(13)
  #assign(paste("perf.a", "1", sep=""),5)
  #perf.a1
  
  thin.data.out.aldex3 <- list() #change name of list here-----------
  data.out.aldex3.u <- list() 
  data.out.aldex3.r <- list() 
  data.out.aldex3.p <- list() 
  
  #for loop
  for (i in 1:nloop){
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin3 <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out.aldex3[[i]] <- thin3
    condsp <- as.vector(thin3$designmat)   # permuted and thinned conditions and data
    datasp <- thin3$mat
    
    #randomized without FP addition PD1
    xrp.aldex3 <- aldex(data, conditions=condsp, gamma = gamma) #uses original dataset but permuted conditions
    data.out.aldex3.r[[i]] <- as.data.frame(xrp.aldex3)
    
    #randomized with FP addition PD1
    xpp.aldex3 <- aldex(datasp, conditions=condsp, gamma = gamma) #uses new dataset with permuted conditions
    data.out.aldex3.p[[i]] <- as.data.frame(xpp.aldex3)
  }
  
  print("done loop")

  #unpermuted PD1
  xup.aldex3 <- aldex(data, conditions=conditions, gamma = gamma)

  return(list(conditions=conditions_p, thin.data=thin.data.out.aldex3, u.data=xup.aldex3, r.data=data.out.aldex3.r, p.data=data.out.aldex3.p))
}
