#function for difference between plots

anal.path <- "../ext_analysis/"
source('code/get_confusion.R')
load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

ald0.conf <- get_confusion(immuno.data_0.aldex2, "ALDEx2", FDR=0.05)
ald3_0.conf <- get_confusion(immuno.data_0.aldex3, "ALDEx3", FDR=0.05)
  
diff.between.1 = ald3_0.conf$diff_coeff[[1]] - ald0.conf$diff_coeff[[1]]
diff.between.2 = ald3_0.conf$diff_coeff[[2]] - ald0.conf$diff_coeff[[2]]
diff.between.3 = ald3_0.conf$diff_coeff[[3]] - ald0.conf$diff_coeff[[3]]
diff.between.4 = ald3_0.conf$diff_coeff[[4]] - ald0.conf$diff_coeff[[4]]
diff.between.5 = ald3_0.conf$diff_coeff[[5]] - ald0.conf$diff_coeff[[5]]
diff.between.6 = ald3_0.conf$diff_coeff[[6]] - ald0.conf$diff_coeff[[6]]
diff.between.7 = ald3_0.conf$diff_coeff[[7]] - ald0.conf$diff_coeff[[7]]
diff.between.8 = ald3_0.conf$diff_coeff[[8]] - ald0.conf$diff_coeff[[8]]
diff.between.9 = ald3_0.conf$diff_coeff[[9]] - ald0.conf$diff_coeff[[9]]
diff.between.10 = ald3_0.conf$diff_coeff[[10]] - ald0.conf$diff_coeff[[10]]
diff.between.11 = ald3_0.conf$diff_coeff[[11]] - ald0.conf$diff_coeff[[11]]
  
ald1.conf <- get_confusion(immuno.data_1.aldex2, "ALDEx2", FDR=0.05)
ald3_1.conf <- get_confusion(immuno.data_1.aldex3, "ALDEx3", FDR=0.05)
  
diff.1between.1 = ald3_1.conf$diff_coeff[[1]] - ald1.conf$diff_coeff[[1]]
diff.1between.2 = ald3_1.conf$diff_coeff[[2]] - ald1.conf$diff_coeff[[2]]
diff.1between.3 = ald3_1.conf$diff_coeff[[3]] - ald1.conf$diff_coeff[[3]]
diff.1between.4 = ald3_1.conf$diff_coeff[[4]] - ald1.conf$diff_coeff[[4]]
diff.1between.5 = ald3_1.conf$diff_coeff[[5]] - ald1.conf$diff_coeff[[5]]
diff.1between.6 = ald3_1.conf$diff_coeff[[6]] - ald1.conf$diff_coeff[[6]]
diff.1between.7 = ald3_1.conf$diff_coeff[[7]] - ald1.conf$diff_coeff[[7]]
diff.1between.8 = ald3_1.conf$diff_coeff[[8]] - ald1.conf$diff_coeff[[8]]
diff.1between.9 = ald3_1.conf$diff_coeff[[9]] - ald1.conf$diff_coeff[[9]]
diff.1between.10 = ald3_1.conf$diff_coeff[[10]] - ald1.conf$diff_coeff[[10]]
diff.1between.11 = ald3_1.conf$diff_coeff[[11]] - ald1.conf$diff_coeff[[11]]
  
ald2.conf <- get_confusion(immuno.data_2.aldex2, "ALDEx2", FDR=0.05)
ald3_2.conf <- get_confusion(immuno.data_2.aldex3, "ALDEx3", FDR=0.05)

diff.2between.1 = ald3_2.conf$diff_coeff[[1]] - ald2.conf$diff_coeff[[1]]
diff.2between.2 = ald3_2.conf$diff_coeff[[2]] - ald2.conf$diff_coeff[[2]]
diff.2between.3 = ald3_2.conf$diff_coeff[[3]] - ald2.conf$diff_coeff[[3]]
diff.2between.4 = ald3_2.conf$diff_coeff[[4]] - ald2.conf$diff_coeff[[4]]
diff.2between.5 = ald3_2.conf$diff_coeff[[5]] - ald2.conf$diff_coeff[[5]]
diff.2between.6 = ald3_2.conf$diff_coeff[[6]] - ald2.conf$diff_coeff[[6]]
diff.2between.7 = ald3_2.conf$diff_coeff[[7]] - ald2.conf$diff_coeff[[7]]
diff.2between.8 = ald3_2.conf$diff_coeff[[8]] - ald2.conf$diff_coeff[[8]]
diff.2between.9 = ald3_2.conf$diff_coeff[[9]] - ald2.conf$diff_coeff[[9]]
diff.2between.10 = ald3_2.conf$diff_coeff[[10]] - ald2.conf$diff_coeff[[10]]
diff.2between.11 = ald3_2.conf$diff_coeff[[11]] - ald2.conf$diff_coeff[[11]]

ald5.conf <- get_confusion(immuno.data_5.aldex2, "ALDEx2", FDR=0.05)
ald3_5.conf <- get_confusion(immuno.data_5.aldex3, "ALDEx3", FDR=0.05)

diff.5between.1 = ald3_5.conf$diff_coeff[[1]] - ald5.conf$diff_coeff[[1]]
diff.5between.2 = ald3_5.conf$diff_coeff[[2]] - ald5.conf$diff_coeff[[2]]
diff.5between.3 = ald3_5.conf$diff_coeff[[3]] - ald5.conf$diff_coeff[[3]]
diff.5between.4 = ald3_5.conf$diff_coeff[[4]] - ald5.conf$diff_coeff[[4]]
diff.5between.5 = ald3_5.conf$diff_coeff[[5]] - ald5.conf$diff_coeff[[5]]
diff.5between.6 = ald3_5.conf$diff_coeff[[6]] - ald5.conf$diff_coeff[[6]]
diff.5between.7 = ald3_5.conf$diff_coeff[[7]] - ald5.conf$diff_coeff[[7]]
diff.5between.8 = ald3_5.conf$diff_coeff[[8]] - ald5.conf$diff_coeff[[8]]
diff.5between.9 = ald3_5.conf$diff_coeff[[9]] - ald5.conf$diff_coeff[[9]]
diff.5between.10 = ald3_5.conf$diff_coeff[[10]] - ald5.conf$diff_coeff[[10]]
diff.5between.11 = ald3_5.conf$diff_coeff[[11]] - ald5.conf$diff_coeff[[11]]
  
  combined0 = list(diff.between.1, diff.between.2, diff.between.3, diff.between.4, diff.between.5, diff.between.6, diff.between.7, diff.between.8, diff.between.9, diff.between.10, diff.between.11)
  combined1 = list(diff.1between.1, diff.1between.2, diff.1between.3, diff.1between.4, diff.1between.5, diff.1between.6, diff.1between.7, diff.1between.8, diff.1between.9, diff.1between.10, diff.1between.11)
  combined2 = list(diff.2between.1, diff.2between.2, diff.2between.3, diff.2between.4, diff.2between.5, diff.2between.6, diff.2between.7, diff.2between.8, diff.2between.9, diff.2between.10, diff.2between.11)
  combined5 = list(diff.5between.1, diff.5between.2, diff.5between.3, diff.5between.4, diff.5between.5, diff.5between.6, diff.5between.7, diff.5between.8, diff.5between.9, diff.5between.10, diff.5between.11)
  
  
  coef <- as.numeric(names(ald0.conf$diff_zero))
  
  
  
  diff.tpr.ald23 <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 10))
  
  for(i in 1:length(coef)){
    temp <- get(paste0("diff.between.", i))
    diff.tpr.ald23[,i] <- abs(temp[,6])
  }
  
  
  #calculations
  coef <- as.numeric(names(ald0.conf$diff_zero))
  TPR_raw <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_raw <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  TPR_0 <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_0 <- matrix(data=NA, nrow=length(ald0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  
  TPR_raw3 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_raw3 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  TPR_03 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  FDR_03 <- matrix(data=NA, nrow=length(ald3_0.conf$diff_zero[[1]]$TP), ncol=length(coef))
  
  for(i in 1:length(coef)){
    TPR_raw[,i] <- ald0.conf$raw_coeff[[i]]$TPR
    FDR_raw[,i] <- ald0.conf$raw_coeff[[i]]$FDR
    TPR_0[,i] <- ald0.conf$raw_zero[[i]]$TPR
    FDR_0[,i] <- ald0.conf$raw_zero[[i]]$FDR
    
    TPR_raw3[,i] <- ald3_0.conf$raw_coeff[[i]]$TPR
    FDR_raw3[,i] <- ald3_0.conf$raw_coeff[[i]]$FDR
    TPR_03[,i] <- ald3_0.conf$raw_zero[[i]]$TPR
    FDR_03[,i] <- ald3_0.conf$raw_zero[[i]]$FDR
  }
  

