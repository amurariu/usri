#line plots

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

#get confusion functions with FDR of 0.1 and 0.05
ald0.conf <- get_confusion(immuno.data_0.aldex2, "ALDEx2", FDR=0.1)
ald1.conf <- get_confusion(aldex_1_var, "ALDEx2", FDR=0.1)
ald2.conf <- get_confusion(aldex_2_var, "ALDEx2", FDR=0.1)
ald5.conf <- get_confusion(aldex_5_var, "ALDEx2", FDR=0.1)
ald3_0.conf <- get_confusion(immuno.data_0.aldex3, "ALDEx3", FDR=0.1)
ald3_1.conf <- get_confusion(aldex3_1_var, "ALDEx3", FDR=0.1)
ald3_2.conf <- get_confusion(aldex3_2_var, "ALDEx3", FDR=0.1)
ald3_5.conf <- get_confusion(aldex3_5_var, "ALDEx3", FDR=0.1)
ald0.conf5 <- get_confusion(aldex_0_var, "ALDEx2", FDR=0.05)
ald1.conf5 <- get_confusion(aldex_1_var, "ALDEx2", FDR=0.05)
ald2.conf5 <- get_confusion(aldex_2_var, "ALDEx2", FDR=0.05)
ald5.conf5 <- get_confusion(aldex_5_var, "ALDEx2", FDR=0.05)
ald3_0.conf5 <- get_confusion(aldex3_0_var, "ALDEx3", FDR=0.05)
ald3_1.conf5 <- get_confusion(aldex3_1_var, "ALDEx3", FDR=0.05)
ald3_2.conf5 <- get_confusion(aldex3_2_var, "ALDEx3", FDR=0.05)
ald3_5.conf5 <- get_confusion(aldex3_5_var, "ALDEx3", FDR=0.05)

#diff betweeen ALDEx2 and ALDEx3

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

combined = list(diff.between.0, diff.between.1, diff.between.2, diff.between.3, diff.between.4, diff.between.5, diff.between.6, diff.between.7, diff.between.8, diff.between.9, diff.between.10)

coef <- as.numeric(names(ald0.conf$diff_zero))

plot(coef, TPR_raw[1,], col=rgb(0,0,0,0.2), ylim=c(0,1), type='l', main = 'ALDEx2')

start <- 1
end <- 100
diff.tpr.ald23 <- as.data.frame(matrix(data = NA, nrow = 1100, ncol = 2))
colnames(diff.tpr.ald23) <- c("coeff", "TPR.difference")

for(i in 1:length(coef)){
  temp <- get(paste0("diff.between.", i))
  diff.tpr.ald23[c(start:end),1] <- rep(coef[i], 100)
  diff.tpr.ald23[c(start:end),2] <- abs(temp[,6])
  start <- start + 100 
  end <- end + 100
}

par(mfrow=c(1,2))

## diff between plot that works

diff.tpr.ald23 <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 10))

for(i in 1:length(coef)){
  temp <- get(paste0("diff.between.", i))
  diff.tpr.ald23[,i] <- abs(temp[,6])
}



par(mfrow=c(1,2))

#ALDEx3 panel
plot(coef, TPR_raw3[1,], col=rgb(0,0.2,0.4,0.2), ylim=c(0,1), type='l', main = 'ALDEx3')
points(coef, FDR_raw3[1,], col=rgb(0,0.6,1,0.2), type='l')
points(coef, TPR_03[1,], col=rgb(1,0.2,0.4,0.2), type='l')
points(coef, FDR_03[1,], col=rgb(0.8,0.2,0.6,0.2), type='l')

for(i in 2:nrow(TPR_raw)){
  points(coef, TPR_raw3[i,], col=rgb(0,0.2,0.4,0.2), type='l')
  points(coef, FDR_raw3[i,], col=rgb(0,0.6,1,0.2), type='l')
  points(coef, TPR_03[i,], col=rgb(1,0.2,0.4,0.2), type='l')
  points(coef, FDR_03[i,], col=rgb(0.8,0.2,0.6,0.2), type='l')
}

#diff.between panel
plot(coef, diff.tpr.ald23[1,], ylim = c(0,0.04), col=rgb(0,0,0,0.2), type='l', main = 'difference between')

for(i in 1:nrow(diff.tpr.ald23)){points(coef, diff.tpr.ald23[i,], ylim = c(0,0.04), col=rgb(0,0,0,0.2), type='l')}


#plots
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
 
par(mfrow=c(1,2))

plot(coef, TPR_raw[1,], col=rgb(0,0,0,0.2), ylim=c(0,1), type='l', main = 'ALDEx2')
points(coef, FDR_raw[1,], col=rgb(0,0,1,0.2), type='l')
points(coef, TPR_0[1,], col=rgb(1,0,0,0.2), type='l')
points(coef, FDR_0[1,], col=rgb(1,0,1,0.2), type='l')

for(i in 2:nrow(TPR_raw)){
points(coef, TPR_raw[i,], col=rgb(0,0,0,0.2), type='l')
points(coef, FDR_raw[i,], col=rgb(0,0,1,0.2), type='l')
points(coef, TPR_0[i,], col=rgb(1,0,0,0.2), type='l')
points(coef, FDR_0[i,], col=rgb(1,0,1,0.2), type='l')
}

plot(coef, TPR_raw3[1,], col=rgb(0,0.2,0.4,0.2), ylim=c(0,1), type='l', main = 'ALDEx3')
points(coef, FDR_raw3[1,], col=rgb(0,0.6,1,0.2), type='l')
points(coef, TPR_03[1,], col=rgb(1,0.2,0.4,0.2), type='l')
points(coef, FDR_03[1,], col=rgb(0.8,0.2,0.6,0.2), type='l')

for(i in 2:nrow(TPR_raw)){
points(coef, TPR_raw3[i,], col=rgb(0,0.2,0.4,0.2), type='l')
points(coef, FDR_raw3[i,], col=rgb(0,0.6,1,0.2), type='l')
points(coef, TPR_03[i,], col=rgb(1,0.2,0.4,0.2), type='l')
points(coef, FDR_03[i,], col=rgb(0.8,0.2,0.6,0.2), type='l')
}