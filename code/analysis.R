anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

source('code/immuno_TPFP.R')

#load files
#immuno dataset
load(paste(anal.path,"immuno.data.deseq.Rda", sep="")) 
load(paste(anal.path,"immuno.data.edger.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_0.out.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.out.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.out.Rda", sep=""))

#inputting parameters for function
immuno.analysis <- analysis.fun(type, , nloop=100) 
save(immuno.analysis, file="./analysis/immuno.data.analysis.Rda") #should it be an Rda file or another file type

  