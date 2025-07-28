source('code/summary.data.fun.R')

load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

sum.fun(aldex2_0 = immuno.data_0.aldex2, aldex2_1 = immuno.data_1.aldex2, aldex2_2 = immuno.data_2.aldex2, aldex2_5 = immuno.data_5.aldex2, aldex3_0 = immuno.data_0.aldex3, aldex3_1 = immuno.data_1.aldex3, aldex3_2 = immuno.data_2.aldex3, aldex3_5 = immuno.data_5.aldex3)

load(paste(anal.path,"lihc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_5.Rda", sep=""))

sum.fun(aldex2_0 = lihc.data_0.aldex2, aldex2_1 = lihc.data_1.aldex2, aldex2_2 = lihc.data_2.aldex2, aldex2_5 = lihc.data_5.aldex2, aldex3_0 = lihc.data_0.aldex3, aldex3_1 = lihc.data_1.aldex3, aldex3_2 = lihc.data_2.aldex3, aldex3_5 = lihc.data_5.aldex3)