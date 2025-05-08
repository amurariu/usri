library(edgeR,quietly=T) 
library(DESeq2,quietly=T)
library(seqgendiff, quietly=T)
library(ALDEx2, warn.conflicts=F)

if(file.exists("https://raw.githubusercontent.com/amurariu/usri/main/analysis/brcadata")){
  load("https://raw.githubusercontent.com/amurariu/usri/main/analysis/brcadata") # file is data.out, list
} else {
  
  raw_counts<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
  con<-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"
    
  brca <- read.table(file=raw_counts, header=T, row.names=1, sep='\t')
  brca.conds <- as.vector(unlist(read.table(file=con, sep='\t')))
  conditions <- data.frame(brca.conds)
  
  ###
  # edgeR functions
  
  y <- DGEList(counts=brca, group=factor(brca.conds))
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  
  # make the filtered base dataset
  brca.data <- y$counts
  brca.data.out <- list()
  
  
  for (i in 1:10){
    
    #unpermuted DESeq2
    dds.u  <- DESeqDataSetFromMatrix(countData = brca.data,
                                     colData = conditions,
                                     design = ~ brca.conds)
    dds.u <- DESeq(dds.u)
    res.u <- results(dds.u)
    
    #permuted + FP addition DESeq2
    thin.brca <- thin_2group(brca.data, prop_null=0.95, alpha=0,
                             signal_fun = stats::rnorm, signal_params = list(mean = 0, sd = 2))
    
    conds <- as.vector(thin.brca$designmat)
    dds.th  <- DESeqDataSetFromMatrix(countData = thin.brca$mat,
                                      colData = data.frame(conds),
                                      design = ~ conds)
    dds.th <- DESeq(dds.th)
    res.th <- results(dds.th)
    
    #unpermuted edgeR
    
    group<-factor(brca.conds)
    design <- model.matrix(~group)
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=2)
    edg.u<-topTags(qlf, n=20478, adjust.method = "BH", sort.by = "none", p.value = 1)
    
    #perm edgeR
    #for me add code here--------------------------
    
    #unperm aldex2
    x <- aldex(brca.data, conditions=brca.conds, gamma=1e-3)
    x.2 <- aldex(brca.data, conditions=brca.conds, gamma=0.2)
    x.5 <- aldex(brca.data, conditions=brca.conds, gamma=0.5)
    
    #perm aldex2
    xp <- aldex(thin.brca$mat, conditions=as.vector(thin.brca$designmat), gamma=1e-3)
    x.2p <- aldex(thin.brca$mat, conditions=as.vector(thin.brca$designmat), gamma=0.2)
    x.5p <- aldex(thin.brca$mat, conditions=as.vector(thin.brca$designmat), gamma=0.5)
    
    data.iter <- list(desu=res.u, desp=res.th, edgeru=edg.u, edgerp=edg.p, ald0u=x, ald2u=x.2, ald5u=x.5, ald0p=xp, ald2p=x.2p, ald5p=x.5p)
    brca.data.out[[i]] <- data.iter
  }
  save(brca.data.out, file="analysis/thin_sim_brca.out.Rda")
}
    
    
  
  
  
  
  
