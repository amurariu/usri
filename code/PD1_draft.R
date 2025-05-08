if(file.exists("https://raw.githubusercontent.com/amurariu/usri/main/analysis/pd1data")){
  load("https://raw.githubusercontent.com/amurariu/usri/main/analysis/pd1data") # file is data.out, list
} else {
  library(ALDEx2, warn.conflicts=F)
  library(seqgendiff, warn.conflicts=F)
  library(edgeR, warn.conflicts=F)
  library(DESeq2, warn.conflicts=F)
  
  raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
  meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
  
  immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
  m <- read.table(file=meta, header=F, row.names=1, sep='\t')
  
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
  res.u <- results(dds.u)
  
  #randomized + TP addition for DESeq2
  
  # do 10 replicates and keep outputs
    # this adds rnorm noise to 5% of the transcripts
    # setting alpha=1 gives no difference between if features are
    # approximately gaussian
    thin.immuno <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                               signal_fun = stats::rnorm, 
                               signal_params = list(mean = 0, sd = 2))
    
    condsp <- as.vector(thin.immuno$designmat)
    
   # setClassUnion("ExpData", c("matrix", "SummarizedExperiment")) #added due to error message being shown for DESeq2, sometimes works and sometimes doesn't?
    dds.th  <- DESeqDataSetFromMatrix(countData = thin.immuno$mat,
                                      colData = data.frame(condsp),
                                      design = ~ condsp)
    dds.th <- DESeq(dds.th)
    res.th <- results(dds.th)
    
    
#unpermuted edgeR
  group<-factor(conditions)
  design <- model.matrix(~group)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  edgeR.res.u<-topTags(qlf, n=nrow(thin.immuno$mat), adjust.method = "BH", sort.by = "none", p.value = 1)

# plot(res.u$padj, edgeR.res.u[[1]]$FDR, log='xy')
  
#randomized + TP addition edgeR - not working
      group_e <- factor(condsp)
      design_e <- model.matrix(~group_e)
      # need to pull from the right slot in thin.immuno
      fit_e <- glmQLFit(thin.immuno$mat,design_e) #negative counts not allowed message
      qlf_e <- glmQLFTest(fit_e,coef=2)
      edgeR.res.p<-topTags(qlf_e, n=nrow(thin.immuno$mat), adjust.method = "BH", sort.by = "none", p.value = 1)
      
      # check to see that things are reasonable
      
     # open circles FP, blue circles TP
     # coef_T <- which(abs(thin.immuno$coefmat) > 0)
     #plot(res.th$padj, edgeR.res.p[[1]]$FDR, log='xy', xlim=c(1e-10,1), ylim=c(1e-10,1))
      #points(res.th$padj[coef_T], edgeR.res.p[[1]]$FDR[coef_T], col='blue', pch=19, cex=0.3)
       
    data.iter <- list(desu=res.u, desp=res.th, edgu=edgeR.res.u, edgp=edgeR.res.p)
    data.out[[i]] <- data.iter
  }
  save(data.out, file="https://raw.githubusercontent.com/amurariu/usri/main/analysis/pd1data")
}

#will create for loop such that it repeats the whole thing and adds it into one file instead of 4 loops










