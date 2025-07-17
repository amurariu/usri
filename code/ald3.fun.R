# data is the raw counts, i.e., immuno from above
# conditions is conditions_p from above
# name is the name of the output file and must be in quotes
# nloops is the number of test loops
ald3.fun <- function(data, conds, gamma, nloop, prop_null, mean){
  
  thin.data.out.aldex3 <- list() 
  aldex.summary.r <- list()
  aldex.summary.t <- list()
  
  #for loop
  for (i in 1:nloop){
    seed = 20 + i
    set.seed(seed)
    
    message(paste0("\nLoop iteration: ", i))
    
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    #generate thin_2group for each dataset as well as labelling for conditions and new dataset
    message('Thinning data...')
    thin3 <- thin_2group(data, prop_null= prop_null, alpha=0,
                        signal_fun = stats::rnorm, 
                        signal_params = list(mean = mean, sd = 2))
    thin.data.out.aldex3[[i]] <- thin3
    conds_th <- as.vector(thin3$designmat)   #thinned conditions
    data_th <- thin3$mat #thinned data
    
    
    # setting conditions for thinned data
    gp1 <- which(conds_th == levels(factor(conds_th))[1])
    gp2 <- which(conds_th == levels(factor(conds_th))[2])
    conditions_t <- conds_th
    conditions_t[gp1] <- 0
    conditions_t[gp2] <- 1
    conditions_t <- as.numeric(conditions_t)
  
    Xt <- formula(~conditions_t)
    datat <- data.frame(conditions_t=conditions_t)
    nsample <- 256
    
    # randomized conditions only (.r)
    message("Running ALDEx3 (scale = ", gamma, ") on original data with randomised groups...")
    data.out.r <- aldex(data, Xt, data=datat, nsample=nsample, scale=clr.sm, gamma=gamma)
    sum.imm.r <- summary.aldex(data.out.r)
    aldex.summary.r[[i]] <- sum.imm.r

    # LFC column is estimate column 3
    # padj column is p.val.adj column 5
    
    
    # randomized and thinned (.t)
    message("Running ALDEx3 (scale = ", gamma, ") on thinned data with randomised groups...")
    data.out.t <- aldex(data_th, Xt, data=datat, nsample=nsample, scale=clr.sm, gamma=gamma)
    sum.imm.t <- summary.aldex(data.out.t)
    aldex.summary.t[[i]] <- sum.imm.t
    
    }
  
  message("\nAll ", i, " iterations of loop finished\n")
  
  #unpermuted conditions (.u)
  # change conditions to binary from named conditions
  gp1 <- which(conds == levels(factor(conds))[1])
  gp2 <- which(conds == levels(factor(conds))[2])
  conditions <- conds
  conditions[gp1] <- 0
  conditions[gp2] <- 1
  conditions <- as.numeric(conditions)
  
  X <- formula(~conditions)
  dataf <- data.frame(conditions=conditions)
  nsample <- 256
  
  message("Running ALDEx3 (scale = ", gamma, ") on original data with non-randomised groups...")
  set.seed(2025)
  data.out.aldex3.u <- aldex(data, X, data=dataf, nsample=nsample, scale=clr.sm, gamma=gamma)
  sum.imm.u <- summary.aldex(data.out.aldex3.u)
  
  return(list(conditions=conditions, thin.data=thin.data.out.aldex3, u.data = sum.imm.u, r.data = aldex.summary.r, t.data = aldex.summary.t))

  }