setwd("/Users/andreeamurariu/Documents/github/usri/")

immuno<-read.table(file="data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv", header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file="data/imm_metadata.txt", header=F, row.names=1, sep='\t')

##stopped editing here, start here

if(file.exists("analysis/thin_sim_data.out.Rda")){
  load("analysis/thin_sim_data.out.Rda") # file is data.out, list
} else {
  library(ALDEx2, warn.conflicts=F)
  library(seqgendiff, warn.conflicts=F)
  library(edgeR, warn.conflicts=F)
  library(DESeq2, warn.conflicts=F)
  
  ###
  # edgeR functions
  conditions <- rep("Pre", 109)
  conditions[grep("_On",m)] <- "On"
  immuno.conds <- data.frame(conditions)
  
  y <- DGEList(counts=immuno, group=factor(conditions))
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  
  # make the filtered base dataset, filters out genes with no expr 
  immuno.data <- y$counts
  
  data.out <- list()
  
  # do 10 replicates and keep outputs
  for(i in 1:10){
    # this adds rnorm noise to 5% of the transcripts
    # setting alpha=1 gives no difference between if features are
    # approximately gaussian
    thin.immuno <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                               signal_fun = stats::rnorm, signal_params = list(mean = 0, sd = 2))
    
    x <- aldex(thin.immuno$mat, conditions=as.vector(thin.immuno$designmat), gamma=1e-3)
    x.2 <- aldex(thin.immuno$mat, conditions=as.vector(thin.immuno$designmat), gamma=0.2)
    x.5 <- aldex(thin.immuno$mat, conditions=as.vector(thin.immuno$designmat), gamma=0.5)
    
    conds <- as.vector(thin.immuno$designmat)
    
    #DESeq2 functions
    dds.th  <- DESeqDataSetFromMatrix(countData = thin.immuno$mat,
                                      colData = data.frame(conds),
                                      design = ~ conds)
    dds.th <- DESeq(dds.th)
    res.th <- results(dds.th)
    data.iter <- list(coef=thin.immuno$coefmat, ald0=x, ald2=x.2, ald5=x.5, des=res.th)
    data.out[[i]] <- data.iter
  }
  save(data.out, file="../analysis/thin_sim_data_draft.out.Rda")
}