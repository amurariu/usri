raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'

immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')

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
  
  # make the filtered base dataset 
  immuno.data <- y$counts
  data.out <- list()
  
  #unpermuted DESeq2
  for (i in 1:10){
  dds.u  <- DESeqDataSetFromMatrix(countData = immuno.data,
                                    colData = immuno.conds,
                                    design = ~ conditions)
  dds.u <- DESeq(dds.u)
  res.u <- results(dds.u)}
  
  #randomized + TP addition for DESeq2
  # do 10 replicates and keep outputs
  for(i in 1:10){
    # this adds rnorm noise to 5% of the transcripts
    # setting alpha=1 gives no difference between if features are
    # approximately gaussian
    thin.immuno <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                               signal_fun = stats::rnorm, signal_params = list(mean = 0, sd = 2))
    
    condsp <- as.vector(thin.immuno$designmat)
    
    setClassUnion("ExpData", c("matrix", "SummarizedExperiment")) #added due to error message being shown for DESeq2
    dds.th  <- DESeqDataSetFromMatrix(countData = thin.immuno$mat,
                                      colData = data.frame(condsp),
                                      design = ~ condsp)
    dds.th <- DESeq(dds.th)
    res.th <- results(dds.th)}
    
    
#unpermuted edgeR
  for (i in 1:10){
  group<-factor(conditions)
  design <- model.matrix(~group)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  edgeR.res.u<-topTags(qlf, n=20478, adjust.method = "BH", sort.by = "none", p.value = 1)}

  
#randomized + TP addition edgeR  
  for (i in 1:10){
    thin.immuno_e <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                               signal_fun = stats::rnorm, signal_params = list(mean = 0, sd = 2)) #confirm if want to use the same thin.immuno group or no
    cond_e <- as.vector(thin.immuno_e$designmat)
    
      group_e<-factor(cond_e)
      design_e <- model.matrix(~group_e)
      fit <- glmQLFit(thin.immuno_e,design_e) #negative counts not allowed message
      
      qlf <- glmQLFTest(fit,coef=2)
      edgeR.res.p<-topTags(qlf, n=20478, adjust.method = "BH", sort.by = "none", p.value = 1)}
      
    #data.iter <- list(coef=thin.immuno$coefmat, ald0=x, ald2=x.2, ald5=x.5, des=res.th)
    #data.out[[i]] <- data.iter
  }
  #save(data.out, file="../analysis/thin_sim_data_draft.out.Rda")
}

#create for loop such that it repeats the whole thing and adds it into one file instead of 4 loops










