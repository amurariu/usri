anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

source('code/immuno_TPFP.R')

#load files
#immuno dataset
load(paste(anal.path,"immuno.data.deseq.Rda", sep="")) 
load(paste(anal.path,"immuno.data.edger.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"immuno.data.limma.Rda", sep=""))

#brca dataset
load(paste(anal.path,"brca.data.deseq.Rda", sep="")) 
load(paste(anal.path,"brca.data.edger.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"brca.data.limma.Rda", sep=""))

#kirc dataset
load(paste(anal.path,"kirc.data.deseq.Rda", sep="")) 
load(paste(anal.path,"kirc.data.edger.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"kirc.data.limma.Rda", sep=""))

#lihc dataset
load(paste(anal.path,"lihc.data.deseq.Rda", sep="")) 
load(paste(anal.path,"lihc.data.edger.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"lihc.data.limma.Rda", sep=""))

# Inputs: "immuno.data_0.aldex2" "immuno.data_2.aldex2" "immuno.data_5.aldex2" "immuno.data.DESeq" "immuno.data.edgeR", "immuno.data.limma"
# Types: "DESeq", "edgeR", "aldex0.we", "aldex0.wi", "aldex0.2.we", "aldex0.2.wi", aldex0.5.we", ""aldex0.5.wi", "limma"

#inputting parameters for function
immuno.deseq.analysis <- analysis.fun(input = immuno.data.DESeq, type = "DESeq", nloop=100)
save(immuno.deseq.analysis, file="./analysis/immuno.deseq.analysis.Rda")

##would we prefer to do aldex wilcoxon or welch's or both?

brca.deseq.analysis <- analysis.fun(input = brca.data.DESeq, type = "DESeq", nloop=100)
save(brca.deseq.analysis, file="./analysis/brca.deseq.analysis.Rda")




  
