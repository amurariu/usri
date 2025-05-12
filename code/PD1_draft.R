library(ALDEx2, warn.conflicts=F)
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

#file.exists name should be changed once code is finalized
##can remove if file exists part
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
    data.o[[i]] <- data.iter
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
  data.out <- list(unpermuted, data.o)
  
  save(data.out, file="./Documents/github/usri/analysis/test.Rda")
}


### analysis draft start from here

analysis.wi <- matrix(data=NA, nrow=300, ncol=7)
analysis.wi <- as.data.frame(analysis.wi)
colnames(analysis.wi) <- c("coeff","iter", "met","PPV","FDR","SEN","SPE")
met = c("ald", "ald2", "ald5", "des", "des5")
analysis.wi[,3] <- rep(met,60)

row=1
for(coeff in c(0.01,0.1,0.2,0.5,0.75,1)){
  for(i in 1:10){
    model <- which(abs(data.out[[i]]$coef) > coeff)
    null.model <- which(abs(data.out[[i]]$coef) < coeff)
    
    TP.wi.ald  <- intersect(which(data.out[[i]]$ald0$wi.eBH < 0.05), model)
    TP.wi.ald.5 <- intersect(which(data.out[[i]]$ald5$wi.eBH < 0.05), model)
    TP.wi.ald.2 <- intersect(which(data.out[[i]]$ald2$wi.eBH < 0.05), model)
    TP.wi.des <- intersect(which(data.out[[i]]$des$padj < 0.05), model)
    TP.wi.des5 <- intersect(which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) >0.5), model)
    
    FN.wi.ald <- setdiff(model, which(data.out[[i]]$ald0$wi.eBH < 0.05))
    FN.wi.ald.5 <- setdiff(model, which(data.out[[i]]$ald5$wi.eBH < 0.05))
    FN.wi.ald.2 <- setdiff(model, which(data.out[[i]]$ald2$wi.eBH < 0.05))
    FN.wi.des <- setdiff(model, which(data.out[[i]]$des$padj < 0.05))
    FN.wi.des5 <- setdiff(model, which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5))
    
    FP.wi.ald <- intersect(which(data.out[[i]]$ald0$wi.eBH < 0.05), null.model)
    FP.wi.ald.2 <- intersect(which(data.out[[i]]$ald2$wi.eBH < 0.05), null.model)
    FP.wi.ald.5 <- intersect(which(data.out[[i]]$ald5$wi.eBH < 0.05), null.model)
    FP.wi.des <- intersect(which(data.out[[i]]$des$padj < 0.05), null.model)
    FP.wi.des5 <- intersect(which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)
    
    TN.wi.ald <- intersect(which(data.out[[i]]$ald0$wi.eBH >= 0.05), null.model)
    TN.wi.ald.2 <- intersect(which(data.out[[i]]$ald2$wi.eBH >= 0.05), null.model)
    TN.wi.ald.5 <- intersect(which(data.out[[i]]$ald5$wi.eBH >= 0.05), null.model)
    TN.wi.des <- intersect(which(data.out[[i]]$des$padj >= 0.05), null.model)
    TN.wi.des5 <- intersect(which(data.out[[i]]$des$padj >= 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)
    
    PPV.wi.ald <- length(TP.wi.ald)/sum(length(TP.wi.ald),length(FP.wi.ald))
    PPV.wi.ald.2 <- length(TP.wi.ald.2)/sum(length(TP.wi.ald.2),length(FP.wi.ald.2))
    PPV.wi.ald.5 <- length(TP.wi.ald.5)/sum(length(TP.wi.ald.5),length(FP.wi.ald.5))
    PPV.wi.des <- length(TP.wi.des)/sum(length(TP.wi.des),length(FP.wi.des))
    PPV.wi.des5 <- length(TP.wi.des5)/sum(length(TP.wi.des5),length(FP.wi.des5))
    
    FDR.wi.ald <- length(FP.wi.ald)/sum(length(TP.wi.ald),length(FP.wi.ald))
    FDR.wi.ald.2 <- length(FP.wi.ald.2)/sum(length(TP.wi.ald.2),length(FP.wi.ald.2))
    FDR.wi.ald.5 <- length(FP.wi.ald.5)/sum(length(TP.wi.ald.5),length(FP.wi.ald.5))
    FDR.wi.des <- length(FP.wi.des)/sum(length(TP.wi.des),length(FP.wi.des))
    FDR.wi.des5 <- length(FP.wi.des5)/sum(length(TP.wi.des5),length(FP.wi.des5))
    
    SEN.wi.ald <- length(TP.wi.ald)/(length(TP.wi.ald) + length(FN.wi.ald))
    SEN.wi.ald.2 <- length(TP.wi.ald.2)/(length(TP.wi.ald.2) + length(FN.wi.ald.2))
    SEN.wi.ald.5 <- length(TP.wi.ald.5)/(length(TP.wi.ald.5) + length(FN.wi.ald.5))
    SEN.wi.des <- length(TP.wi.des)/(length(TP.wi.des) + length(FN.wi.des))
    SEN.wi.des5 <- length(TP.wi.des5)/(length(TP.wi.des5) + length(FN.wi.des5))
    
    SPE.wi.ald <- length(TN.wi.ald)/(length(TN.wi.ald) + length(FP.wi.ald))
    SPE.wi.ald.2 <- length(TN.wi.ald.2)/(length(TN.wi.ald.2) + length(FP.wi.ald.2))
    SPE.wi.ald.5 <- length(TN.wi.ald.5)/(length(TN.wi.ald.5) + length(FP.wi.ald.5))
    SPE.wi.des <- length(TN.wi.des)/(length(TN.wi.des) + length(FP.wi.des))
    SPE.wi.des5 <- length(TN.wi.des5)/(length(TN.wi.des5) + length(FP.wi.des5))
    
    # even butt uglier ...
    for (j in 1:5){
      analysis.wi[row,1] <- coeff
      analysis.wi[row,2] <- i
      if (j==1){
        analysis.wi[row,4] <- PPV.wi.ald
        analysis.wi[row,5] <- FDR.wi.ald
        analysis.wi[row,6] <- SEN.wi.ald
        analysis.wi[row,7] <- SPE.wi.ald
      } else if (j==2) {
        analysis.wi[row,4] <- PPV.wi.ald.2
        analysis.wi[row,5] <- FDR.wi.ald.2
        analysis.wi[row,6] <- SEN.wi.ald.2
        analysis.wi[row,7] <- SPE.wi.ald.2
      } else if (j==3){
        analysis.wi[row,4] <- PPV.wi.ald.5
        analysis.wi[row,5] <- FDR.wi.ald.5
        analysis.wi[row,6] <- SEN.wi.ald.5
        analysis.wi[row,7] <- SPE.wi.ald.5
      } else if (j==4){
        analysis.wi[row,4] <- PPV.wi.des
        analysis.wi[row,5] <- FDR.wi.des
        analysis.wi[row,6] <- SEN.wi.des
        analysis.wi[row,7] <- SPE.wi.des
      } else if (j==5){
        analysis.wi[row,4] <- PPV.wi.des5
        analysis.wi[row,5] <- FDR.wi.des5
        analysis.wi[row,6] <- SEN.wi.des5
        analysis.wi[row,7] <- SPE.wi.des5
      }
      row=row+1
    }
  }
}

means.wi <- matrix(data=NA, nrow=30, ncol=6)
means.wi <- as.data.frame(means.wi)
colnames(means.wi) <- c("coef", "method", "PPV","FDR","SEN","SPE")

for(i in 1:5){
  means.wi[i,1:2] <- analysis.wi[seq(from=i, to=50, by=5),c(1,3)][1,]
  means.wi[i,3:6] <- colMeans(analysis.wi[seq(from=i, to=50, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+5,1:2] <- analysis.wi[seq(from=i+50, to=100, by=5),c(1,3)][1,]
  means.wi[i+5,3:6] <- colMeans(analysis.wi[seq(from=i+50, to=100, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+10,1:2] <- analysis.wi[seq(from=i+100, to=150, by=5),c(1,3)][1,]
  means.wi[i+10,3:6] <- colMeans(analysis.wi[seq(from=i+100, to=150, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+15,1:2] <- analysis.wi[seq(from=i+150, to=200, by=5),c(1,3)][1,]
  means.wi[i+15,3:6] <- colMeans(analysis.wi[seq(from=i+150, to=200, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+20,1:2] <- analysis.wi[seq(from=i+201, to=250, by=5),c(1,3)][1,]
  means.wi[i+20,3:6] <- colMeans(analysis.wi[seq(from=i+201, to=250, by=5),4:7])
}
for(i in 1:5){
  means.wi[i+25,1:2] <- analysis.wi[seq(from=i+251, to=300, by=5),c(1,3)][1,]
  means.wi[i+25,3:6] <- colMeans(analysis.wi[seq(from=i+251, to=300, by=5),4:7])
}

#### we.eBH
analysis.out <- matrix(data=NA, nrow=300, ncol=7)
analysis.out <- as.data.frame(analysis.out)
colnames(analysis.out) <- c("coeff","iter", "met","PPV","FDR","SEN","SPE")
met = c("ald", "ald2", "ald5", "des", "des5")
analysis.out[,3] <- rep(met,60)

row=1
for(coeff in c(0.01,0.1,0.2,0.5,0.75,1)){
  for(i in 1:10){
    model <- which(abs(data.out[[i]]$coef) > coeff)
    null.model <- which(abs(data.out[[i]]$coef) < coeff)
    
    TP.ald  <- intersect(which(data.out[[i]]$ald0$we.eBH < 0.05), model)
    TP.ald.5 <- intersect(which(data.out[[i]]$ald5$we.eBH < 0.05), model)
    TP.ald.2 <- intersect(which(data.out[[i]]$ald2$we.eBH < 0.05), model)
    TP.des <- intersect(which(data.out[[i]]$des$padj < 0.05), model)
    TP.des5 <- intersect(which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) >0.5), model)
    
    FN.ald <- setdiff(model, which(data.out[[i]]$ald0$we.eBH < 0.05))
    FN.ald.5 <- setdiff(model, which(data.out[[i]]$ald5$we.eBH < 0.05))
    FN.ald.2 <- setdiff(model, which(data.out[[i]]$ald2$we.eBH < 0.05))
    FN.des <- setdiff(model, which(data.out[[i]]$des$padj < 0.05))
    FN.des5 <- setdiff(model, which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5))
    
    FP.ald <- intersect(which(data.out[[i]]$ald0$we.eBH < 0.05), null.model)
    FP.ald.2 <- intersect(which(data.out[[i]]$ald2$we.eBH < 0.05), null.model)
    FP.ald.5 <- intersect(which(data.out[[i]]$ald5$we.eBH < 0.05), null.model)
    FP.des <- intersect(which(data.out[[i]]$des$padj < 0.05), null.model)
    FP.des5 <- intersect(which(data.out[[i]]$des$padj < 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)
    
    TN.ald <- intersect(which(data.out[[i]]$ald0$we.eBH >= 0.05), null.model)
    TN.ald.2 <- intersect(which(data.out[[i]]$ald2$we.eBH >= 0.05), null.model)
    TN.ald.5 <- intersect(which(data.out[[i]]$ald5$we.eBH >= 0.05), null.model)
    TN.des <- intersect(which(data.out[[i]]$des$padj >= 0.05), null.model)
    TN.des5 <- intersect(which(data.out[[i]]$des$padj >= 0.05 & abs(data.out[[i]]$des$log2FoldChange) > 0.5), null.model)
    
    PPV.ald <- length(TP.ald)/sum(length(TP.ald),length(FP.ald))
    PPV.ald.2 <- length(TP.ald.2)/sum(length(TP.ald.2),length(FP.ald.2))
    PPV.ald.5 <- length(TP.ald.5)/sum(length(TP.ald.5),length(FP.ald.5))
    PPV.des <- length(TP.des)/sum(length(TP.des),length(FP.des))
    PPV.des5 <- length(TP.des5)/sum(length(TP.des5),length(FP.des5))
    
    FDR.ald <- length(FP.ald)/sum(length(TP.ald),length(FP.ald))
    FDR.ald.2 <- length(FP.ald.2)/sum(length(TP.ald.2),length(FP.ald.2))
    FDR.ald.5 <- length(FP.ald.5)/sum(length(TP.ald.5),length(FP.ald.5))
    FDR.des <- length(FP.des)/sum(length(TP.des),length(FP.des))
    FDR.des5 <- length(FP.des5)/sum(length(TP.des5),length(FP.des5))
    
    SEN.ald <- length(TP.ald)/(length(TP.ald) + length(FN.ald))
    SEN.ald.2 <- length(TP.ald.2)/(length(TP.ald.2) + length(FN.ald.2))
    SEN.ald.5 <- length(TP.ald.5)/(length(TP.ald.5) + length(FN.ald.5))
    SEN.des <- length(TP.des)/(length(TP.des) + length(FN.des))
    SEN.des5 <- length(TP.des5)/(length(TP.des5) + length(FN.des5))
    
    SPE.ald <- length(TN.ald)/(length(TN.ald) + length(FP.ald))
    SPE.ald.2 <- length(TN.ald.2)/(length(TN.ald.2) + length(FP.ald.2))
    SPE.ald.5 <- length(TN.ald.5)/(length(TN.ald.5) + length(FP.ald.5))
    SPE.des <- length(TN.des)/(length(TN.des) + length(FP.des))
    SPE.des5 <- length(TN.des5)/(length(TN.des5) + length(FP.des5))
    

