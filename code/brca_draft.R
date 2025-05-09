library(edgeR,quietly=T) 
library(DESeq2,quietly=T)
library(seqgendiff, quietly=T)
library(ALDEx2, warn.conflicts=F)

if(file.exists(file="./Documents/github/usri/analysis/test2.Rda")){
  load(file="./Documents/github/usri/analysis/test2.Rda") # file is data.out, list
} else {
  
  raw_counts<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
  con<-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"
    
  brca <- read.table(file=raw_counts, header=T, row.names=1, sep='\t')
  brca.conds <- as.vector(unlist(read.table(file=con, sep='\t')))
  conditions <- data.frame(brca.conds)
  
  # edgeR functions
  y <- DGEList(counts=brca, group=factor(brca.conds))
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  
  # make the filtered base dataset
  brca.data <- y$counts
  brca.data.out <- list()
  
  
  for (i in 1:2){
    
    #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
    thin.brca <- thin_2group(brca.data, prop_null=0.95, alpha=0,
                             signal_fun = stats::rnorm, signal_params = list(mean = 0, sd = 2))
    
    # permuted and thinned conditions and data
    conds <- as.vector(thin.brca$designmat)
    datasp <- thin.brca$mat
    
    #randomized but no FP addition DESeq2
    dds.r  <- DESeqDataSetFromMatrix(countData = brca.data,  #uses original data (no TP added)
                                     colData = data.frame(conds), #uses data randomization order from thin.brca
                                     design = ~ conds)
    dds.r <- DESeq(dds.r)
    res.r <- results(dds.r)
    
    #permuted + FP addition DESeq2
    dds.th  <- DESeqDataSetFromMatrix(countData = thin.brca$mat, #data that includes TP
                                      colData = data.frame(conds),
                                      design = ~ conds)
    dds.th <- DESeq(dds.th)
    res.th <- results(dds.th)
    
    #randomized but no FP addition edgeR
    group_e <- factor(conds)
    design_e <- model.matrix(~group_e) #use data randomization from seqgendiff
    fit_ <- glmQLFit(brca.data,design_e) #uses original data (ie. no TP added)
    qlf_ <- glmQLFTest(fit_,coef=2)
    edg.r<-topTags(qlf_, n=nrow(brca.data), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    
    #permuted + FP addition edgeR
    fit_e <- glmQLFit(datasp,design_e) #data with TPs
    qlf_e <- glmQLFTest(fit_e,coef=2)
    edg.p<-topTags(qlf_e, n=nrow(datasp), adjust.method = "BH", sort.by = "none", p.value = 1)
    
    #permuted without TP addition aldex2
    xr <- aldex(brca.data, conditions=conds, gamma=1e-3) #uses original dataset but permuted conditions
    x.2r <- aldex(brca.data, conditions=conds, gamma=0.2)
    x.5r <- aldex(brca.data, conditions=conds, gamma=0.5)
   
    
     #permuted + TP addition aldex2
    xp <- aldex(datasp, conditions=conds, gamma=1e-3)
    x.2p <- aldex(datasp, conditions=conds, gamma=0.2)
    x.5p <- aldex(datasp, conditions=conds, gamma=0.5)
    
    data.iter <- list(desr=res.r, desp=res.th, edgr=edg.r, edgp=edg.p,ald0r=xr, ald2r=x.2r, ald5r=x.5r, ald0p=xp, ald2p=x.2p, ald5p=x.5p )
    brca.data.out[[i]] <- data.iter
  }
  
  #unpermuted DESeq2
  dds.u  <- DESeqDataSetFromMatrix(countData = brca.data,
                                   colData = conditions,
                                   design = ~ brca.conds)
  dds.u <- DESeq(dds.u)
  res.u <- results(dds.u)
  
  #unpermuted edgeR
  group<-factor(brca.conds)
  design <- model.matrix(~group)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  edg.u<-topTags(qlf, n=nrow(brca.data), adjust.method = "BH", sort.by = "none", p.value = 1)
  
  #unpermuted aldex2
  x <- aldex(brca.data, conditions=brca.conds, gamma=1e-3)
  x.2 <- aldex(brca.data, conditions=brca.conds, gamma=0.2)
  x.5 <- aldex(brca.data, conditions=brca.conds, gamma=0.5)

  unpermuted<-list(desu=res.u, edgeru=edg.u, ald0u=x, ald2u=x.2, ald5u=x.5)
  combinedbrca <- list(unpermuted, brca.data.out)
  
  save(combinedbrca, file="./Documents/github/usri/analysis/test2.Rda")
}
    
    
  
  
  
  
  
