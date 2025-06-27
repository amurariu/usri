#testing graph.fun 

source('code/graph.fun.R')

anal.path <- "../ext_analysis/" #move ext_analysis back to github folder

#load files
#immuno dataset
load(paste(anal.path,"immuno.data.deseq.Rda", sep="")) 
load(paste(anal.path,"immuno.data.edger.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"immuno.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"immuno.data.limma.Rda", sep=""))

graph.fun(DESeq_var = immuno.data.DESeq, edgeR_var = immuno.data.edgeR, limma_var = immuno.data.limma, aldex_0_var = immuno.data_0.aldex2, aldex_2_var = immuno.data_2.aldex2, aldex_5_var = immuno.data_5.aldex2)

#brca dataset
load(paste(anal.path,"brca.data.deseq.Rda", sep="")) 
load(paste(anal.path,"brca.data.edger.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"brca.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"brca.data.limma.Rda", sep=""))

graph.fun(DESeq_var = brca.data.DESeq, edgeR_var = brca.data.edgeR, limma_var = brca.data.limma, aldex_0_var = brca.data_0.aldex2, aldex_2_var = brca.data_2.aldex2, aldex_5_var = brca.data_5.aldex2)


#kirc dataset
load(paste(anal.path,"kirc.data.deseq.Rda", sep="")) 
load(paste(anal.path,"kirc.data.edger.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"kirc.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"kirc.data.limma.Rda", sep=""))

graph.fun(DESeq_var = kirc.data.DESeq, edgeR_var = kirc.data.edgeR, limma_var = kirc.data.limma, aldex_0_var = kirc.data_0.aldex2, aldex_2_var = kirc.data_2.aldex2, aldex_5_var = kirc.data_5.aldex2)


#lihc dataset
load(paste(anal.path,"lihc.data.deseq.Rda", sep="")) 
load(paste(anal.path,"lihc.data.edger.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"lihc.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"lihc.data.limma.Rda", sep=""))

graph.fun(DESeq_var = lihc.data.DESeq, edgeR_var = lihc.data.edgeR, limma_var = lihc.data.limma, aldex_0_var = lihc.data_0.aldex2, aldex_2_var = lihc.data_2.aldex2, aldex_5_var = lihc.data_5.aldex2)

#luad dataset
load(paste(anal.path,"luad.data.deseq.Rda", sep="")) 
load(paste(anal.path,"luad.data.edger.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"luad.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"luad.data.limma.Rda", sep=""))

graph.fun(DESeq_var = luad.data.DESeq, edgeR_var = luad.data.edgeR, limma_var = luad.data.limma, aldex_0_var = luad.data_0.aldex2, aldex_2_var = luad.data_2.aldex2, aldex_5_var = luad.data_5.aldex2)


#prad dataset
load(paste(anal.path,"prad.data.deseq.Rda", sep="")) 
load(paste(anal.path,"prad.data.edger.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"prad.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"prad.data.limma.Rda", sep=""))

graph.fun(DESeq_var = prad.data.DESeq, edgeR_var = prad.data.edgeR, limma_var = prad.data.limma, aldex_0_var = prad.data_0.aldex2, aldex_2_var = prad.data_2.aldex2, aldex_5_var = prad.data_5.aldex2)

#thca dataset
load(paste(anal.path,"thca.data.deseq.Rda", sep="")) 
load(paste(anal.path,"thca.data.edger.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_0.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_2.Rda", sep=""))
load(paste(anal.path,"thca.data.aldex2_5.Rda", sep=""))
load (paste(anal.path,"thca.data.limma.Rda", sep=""))

graph.fun(DESeq_var = thca.data.DESeq, edgeR_var = thca.data.edgeR, limma_var = thca.data.limma, aldex_0_var = thca.data_0.aldex2, aldex_2_var = thca.data_2.aldex2, aldex_5_var = thca.data_5.aldex2)

