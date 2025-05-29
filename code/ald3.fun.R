# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald3.fun <- function(Y, data, conditions, nsample, nloop=1, scale, gamma){
  set.seed(4)
  #assign(paste("perf.a", "1", sep=""),5)
  #perf.a1
  #gam <- c(1e-3, 0.2, 0.5) #contains different scale values
  
  thin.data.out.aldex3 <- list() #change name of list here-----------
  data.out.aldex.u3 <- list() 
  data.out.aldex.r3 <- list() 
  data.out.aldex.p3 <- list() 
  
  #for loop
  for (i in 1:nloop){
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin3 <- thin_2group(data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out.aldex3[[i]] <- thin3
    condsp <- as.vector(thin3$designmat) # permuted and thinned conditions and data
    datasp <- thin3$mat
    Y_datasp <- matrix(nrow=nrow(datasp), ncol=ncol(datasp))
    
    #randomized without FP addition PD1
    xrp.aldex3 <- aldex(Y, data, conditions=condsp, nsample=1024, scale = scale, gamma = gamma) #uses original dataset but permuted conditions
    data.out.aldex.r3[[i]] <- xrp.aldex3
    
    #randomized with FP addition PD1
    xpp.aldex3 <- aldex(Y = Y_datasp, datasp, conditions=condsp, nsample=1024, scale = clr.sm, gamma = gamma) #uses new dataset with permuted conditions
    data.out.aldex.p3[[i]] <- xpp.aldex3
  }
  
  print("done loop")
  
  #unpermuted PD1
  xup.aldex <- aldex(Y, data, conditions=conditions, nsample=1024, gamma = gamma)
  
  return(list(conditions=conditions, thin.data=thin.data.out.aldex3, u.data=xup.aldex, r.data=data.out.aldex.r3, p.data=data.out.aldex.p3))
}
