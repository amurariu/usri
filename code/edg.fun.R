# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
edg.fun <- function(data, conditions, nloop=4){
  
  #assign(paste("perf.a", "1", sep=""),5)
  #perf.a1
  conditions_p <- conditions
  conds <- data.frame(conditions_p)
  
  thin.data.out.edger <- list() #change name of list here-----------
  data.out.edgeR.u <- list() 
  data.out.edgeR.r <- list() 
  data.out.edgeR.p <- list() 
  
  #for loop
  for (i in 1:nloop){
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    
    #PD1
    thin <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                               signal_fun = stats::rnorm, 
                               signal_params = list(mean = 0, sd = 2))
    thin.data.out.edger[[i]] <- thin
    condsp <- as.vector(thin$designmat)   # permuted and thinned conditions and data
    datasp <- thin$mat
    
    #edgeR analysis
    #PD1 setup
    group_p <- factor(condsp)
    design_p <- model.matrix(~group_p) #use data randomization from seqgendiff
    
    #randomized without FP addition PD1
    fit_rp <- glmQLFit(immuno.data,design_p) #uses original data (ie. no TP added)
    qlf_rp <- glmQLFTest(fit_rp,coef=2)
    edg.rp<-topTags(qlf_rp, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    resrp.edgeR<-list(resu=edg.rp)
    data.out.edgeR.r[[i]] <- as.data.frame(resrp.edgeR)
    
    #randomized with FP addition PD1
    fit_pp <- glmQLFit(datasp,design_p)
    qlf_pp <- glmQLFTest(fit_pp,coef=2)
    edg.pp<-topTags(qlf_pp, n=nrow(datasp), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    respp.edgeR<-list(resu=edg.pp)
    data.out.edgeR.p[[i]] <- as.data.frame(respp.edgeR)
  }
  print("done loop")
  
  
  #unpermuted PD1
  group_up<-factor(conditions_p)
  design_up <- model.matrix(~group_up)
  fit_up <- glmQLFit(y_pd1,design_up)
  qlf_up <- glmQLFTest(fit_up,coef=2)
  edg.up<-topTags(qlf_up, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1) 
  
  resup.edgeR<-list(resu=edg.up)
  data.out.edgeR.u <- list(resup.edgeR)
  
  return(list(conditions=conditions_p, thin.data=thin.data.out, u.data=data.edg.u, r.data=data.out.edger.r, p.data=data.out.edger.p))
}
