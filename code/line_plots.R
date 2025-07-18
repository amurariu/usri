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

diff.between.0 = unlist(ald3_0.conf) - unlist(ald0.conf)

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