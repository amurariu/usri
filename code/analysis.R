anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

source('code/immuno_TPFP.R')

#load files
#immuno dataset
load(paste(anal.path,"immuno.data.deseq.Rda", sep="")) 
load(paste(anal.path,"immuno.data.edger.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_0.out.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.out.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.out.Rda", sep=""))

# Inputs: "immuno.data_0.aldex2" "immuno.data_2.aldex2" "immuno.data_5.aldex2" "immuno.data.DESeq" "immuno.data.edgeR"
# Types: "DESeq", "edgeR", "aldex0.we", "aldex0.wi", "aldex0.2.we", "aldex0.2.wi", aldex0.5.we", ""aldex0.5.wi"

#inputting parameters for function
immuno.deseq.analysis <- analysis.fun(input = "immuno.data.DESeq", type = "DESeq", nloop=1)
save(immuno.deseq.analysis, file="./analysis/immuno.deseq.analysis.Rda")

  