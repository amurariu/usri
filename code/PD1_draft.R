library(ALDEx2, warn.conflicts=F)
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

#file.exists name should be changed once code is finalized
if(file.exists("https://raw.githubusercontent.com/amurariu/usri/main/analysis/test.Rda")){
  load("https://raw.githubusercontent.com/amurariu/usri/main/analysis/test.Rda") # file is data.out, list
} else {
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

  for (i in 1:2){
  # do 10 replicates and keep outputs
    # this adds rnorm noise to 5% of the transcripts
    # setting alpha=1 gives no difference between if features are
    # approximately gaussian
    thin.immuno <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                               signal_fun = stats::rnorm, 
                               signal_params = list(mean = 0, sd = 2))
    # permuted and thinned conditions and data
    condsp <- as.vector(thin.immuno$designmat)
    datasp <- thin.immuno$mat
    
    
     #randomized DESeq2 with no TP generation
    
     dds.r  <- DESeqDataSetFromMatrix(countData = immuno.data,  #uses original data (no TP added)
                                       colData = data.frame(condsp), #uses data randomization order from thin.immuno
                                       design = ~ condsp)
     dds.r <- DESeq(dds.r)
     res.r <- results(dds.r)
    
    
     #randomized + TP addition for DESeq2
    dds.th  <- DESeqDataSetFromMatrix(countData = datasp,
                                      colData = data.frame(condsp),
                                      design = ~ condsp)
    dds.th <- DESeq(dds.th)
    res.th <- results(dds.th)

    # plot(res.u$padj, edgeR.res.u[[1]]$FDR, log='xy')
  
    #randomized edgeR with no TP generation
    group_e <- factor(condsp)
    design_e <- model.matrix(~group_e) #use data randomization from seqgendiff
    # need to pull from the right slot in thin.immuno
    fit_ <- glmQLFit(immuno.data,design_e) #uses original data (ie. no TP added)
    qlf_ <- glmQLFTest(fit_,coef=2)
    edgeR.res.r<-topTags(qlf_, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    
    #randomized + TP addition edgeR
      fit_e <- glmQLFit(datasp,design_e)
      qlf_e <- glmQLFTest(fit_e,coef=2)
      edgeR.res.p<-topTags(qlf_e, n=nrow(datasp), adjust.method = "BH", sort.by = "none", p.value = 1)
      
     # check to see that things are reasonable
     # open circles FP, blue circles TP
     # coef_T <- which(abs(thin.immuno$coefmat) > 0)
     #plot(res.th$padj, edgeR.res.p[[1]]$FDR, log='xy', xlim=c(1e-10,1), ylim=c(1e-10,1))
      #points(res.th$padj[coef_T], edgeR.res.p[[1]]$FDR[coef_T], col='blue', pch=19, cex=0.3)
      
      #permuted without TP addition aldex2
      xr <- aldex(immuno.data, conditions=condsp, gamma=1e-3) #uses original dataset but permuted conditions
      x.2r <- aldex(immuno.data, conditions=condsp, gamma=0.2)
      x.5r <- aldex(immuno.data, conditions=condsp, gamma=0.5)
      
      
      #permuted + TP addition aldex2
      xp <- aldex(datasp, conditions=condsp, gamma=1e-3)
      x.2p <- aldex(datasp, conditions=condsp, gamma=0.2)
      x.5p <- aldex(datasp, conditions=condsp, gamma=0.5)
       
    data.iter <- list(desp=res.th, desr=res.r, edgp=edgeR.res.p, edgr= edgeR.res.r, ald0r=xr, ald2r=x.2r, ald5r=x.5r, ald0p=xp, ald2p=x.2p, ald5p=x.5p)
    data.out[[i]] <- data.iter
  }
  
  #unpermuted DESeq2
  
  #setClassUnion("ExpData", c("matrix", "SummarizedExperiment")) #added due to error message being shown for DESeq2, sometimes works and sometimes doesn't?
  dds.u  <- DESeqDataSetFromMatrix(countData = immuno.data,
                                   colData = immuno.conds,
                                   design = ~ conditions)
  dds.u <- DESeq(dds.u)
  res.u <- results(dds.u)
  
  #unpermuted edgeR
  group<-factor(conditions)
  design <- model.matrix(~group)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  edgeR.res.u<-topTags(qlf, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1) 
  
  #unpermuted aldex2
  x <- aldex(immuno.data, conditions=immuno.conds, gamma=1e-3)
  x.2 <- aldex(immuno.data, conditions=immuno.conds, gamma=0.2)
  x.5 <- aldex(immuno.data, conditions=immuno.conds, gamma=0.5)
  
  unpermuted<-list(desu=res.u, edgu=edgeR.res.u, ald0u=x, ald2u=x.2, ald5u=x.5)
  combined <- list(unpermuted, data.out)
  
  save(combined, file="./Documents/github/usri/analysis/test.Rda")
}
