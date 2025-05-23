# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald.fun <- function(data, conditions, nloop=4){
  
  #assign(paste("perf.a", "1", sep=""),5)
  #perf.a1
  conditions_p <- conditions
  conds <- data.frame(conditions_p)
  
  thin.data.out.aldex <- list() #change name of list here-----------
  data.out.aldex.u <- list() 
  data.out.aldex.r <- list() 
  data.out.aldex.p <- list() 
  
  #for loop
  for (i in 1:nloop){
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = 0, sd = 2))
    thin.data.out.aldex[[i]] <- thin
    condsp <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    datasp <- thin$mat
    
    #randomized without FP addition PD1
    xrp.aldex <- aldex(immuno.data, conditions=condsp, gamma=1e-3) #uses original dataset but permuted conditions
    resrp.aldex<-list(resu=xrp.aldex)
    data.out.aldex.r[[i]] <- as.data.frame(resrp.aldex)
    
    #randomized with FP addition PD1
    xpp.aldex <- aldex(datasp, conditions=condsp, gamma=1e-3) #uses new dataset with permuted conditions
    respp.aldex<-list(resu=xpp.aldex)
    data.out.aldex.p[[i]] <- as.data.frame(respp.aldex)
  }
  print("done loop")
  
  #unpermuted PD1
  xup.aldex <- aldex(immuno.data, conditions=immuno.conds, gamma=1e-3)
  immuno.data.out.aldex.u <- list(xup.aldex)
  
  return(list(conditions=conditions_p, thin.data=thin.data.out.aldex, u.data=data.out.aldex.u, r.data=data.out.aldex.r, p.data=data.out.aldex.p))
}
