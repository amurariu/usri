anal.path <- "../ext_analysis/"
source('code/summary.data.fun.R')

load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex3_5.Rda", sep=""))

imm.sum <- sum.fun(aldex2_0 = immuno.data_0.aldex2, aldex2_1 = immuno.data_1.aldex2, aldex2_2 = immuno.data_2.aldex2, aldex2_5 = immuno.data_5.aldex2, aldex3_0 = immuno.data_0.aldex3, aldex3_1 = immuno.data_1.aldex3, aldex3_2 = immuno.data_2.aldex3, aldex3_5 = immuno.data_5.aldex3)
save(imm.sum, file="analysis/summarystats/ss.immuno.Rda")


#cannot load all files at once
load(paste(anal.path,"brca.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex3_5.Rda", sep=""))

brca.sum <- sum.fun(aldex2_0 = brca.data_0.aldex2, aldex2_1 = brca.data_1.aldex2, aldex2_2 = brca.data_2.aldex2, aldex2_5 = brca.data_5.aldex2, aldex3_0 = brca.data_0.aldex3, aldex3_1 = brca.data_1.aldex3, aldex3_2 = brca.data_2.aldex3, aldex3_5 = brca.data_5.aldex3)
save(brca.sum, file="analysis/summarystats/ss.brca.Rda")

#load this next
load(paste(anal.path,"kirc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex3_5.Rda", sep=""))

kirc.sum <- sum.fun(aldex2_0 = kirc.data_0.aldex2, aldex2_1 = kirc.data_1.aldex2, aldex2_2 = kirc.data_2.aldex2, aldex2_5 = kirc.data_5.aldex2, aldex3_0 = kirc.data_0.aldex3, aldex3_1 = kirc.data_1.aldex3, aldex3_2 = kirc.data_2.aldex3, aldex3_5 = kirc.data_5.aldex3)
save(kirc.sum, file="analysis/summarystats/ss.kirc.Rda")


load(paste(anal.path,"lihc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex3_5.Rda", sep=""))

lihc.sum <- sum.fun(aldex2_0 = lihc.data_0.aldex2, aldex2_1 = lihc.data_1.aldex2, aldex2_2 = lihc.data_2.aldex2, aldex2_5 = lihc.data_5.aldex2, aldex3_0 = lihc.data_0.aldex3, aldex3_1 = lihc.data_1.aldex3, aldex3_2 = lihc.data_2.aldex3, aldex3_5 = lihc.data_5.aldex3)
save(lihc.sum, file="analysis/summarystats/ss.lihc.Rda")

load(paste(anal.path,"luad.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex3_5.Rda", sep=""))

luad.sum <- sum.fun(aldex2_0 = luad.data_0.aldex2, aldex2_1 = luad.data_1.aldex2, aldex2_2 = luad.data_2.aldex2, aldex2_5 = luad.data_5.aldex2, aldex3_0 = luad.data_0.aldex3, aldex3_1 = luad.data_1.aldex3, aldex3_2 = luad.data_2.aldex3, aldex3_5 = luad.data_5.aldex3)
save(luad.sum, file="analysis/summarystats/ss.luad.Rda")

load(paste(anal.path,"prad.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex3_5.Rda", sep=""))

prad.sum <- sum.fun(aldex2_0 = prad.data_0.aldex2, aldex2_1 = prad.data_1.aldex2, aldex2_2 = prad.data_2.aldex2, aldex2_5 = prad.data_5.aldex2, aldex3_0 = prad.data_0.aldex3, aldex3_1 = prad.data_1.aldex3, aldex3_2 = prad.data_2.aldex3, aldex3_5 = prad.data_5.aldex3)
save(prad.sum, file="analysis/summarystats/ss.prad.Rda")

load(paste(anal.path,"thca.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_1.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_5.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_0.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_1.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_2.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex3_5.Rda", sep=""))

thca.sum <- sum.fun(aldex2_0 = thca.data_0.aldex2, aldex2_1 = thca.data_1.aldex2, aldex2_2 = thca.data_2.aldex2, aldex2_5 = thca.data_5.aldex2, aldex3_0 = thca.data_0.aldex3, aldex3_1 = thca.data_1.aldex3, aldex3_2 = thca.data_2.aldex3, aldex3_5 = thca.data_5.aldex3)
save(thca.sum, file="analysis/summarystats/ss.thca.Rda")



