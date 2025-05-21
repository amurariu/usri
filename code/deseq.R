library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

#download datasets

#data_files = c("immuno", "brca")
#for (i in 1:data_files) {

#immuno/PD1 dataset loading
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p)


#edgeR conditions for initial filtering
#PD1
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset
immuno.data.out.deseq.u <- list() 
immuno.data.out.deseq.r <- list() 
immuno.data.out.deseq.p <- list() 

#for loop
for (i in 1:2){
  #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
  #generate thin_2group for each dataset as well as labelling for conditions and new dataset
  
  #PD1
  thin.immuno <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                             signal_fun = stats::rnorm, 
                             signal_params = list(mean = 0, sd = 2))
  condsp <- as.vector(thin.immuno$designmat)   # permuted and thinned conditions and data
  datasp <- thin.immuno$mat
  
  #DESeq2 analysis
  #randomized without FP addition PD1
  dds.rp.deseq  <- DESeqDataSetFromMatrix(countData = immuno.data,  #uses original data (no TP added)
                                   colData = data.frame(condsp), #uses data randomization order from thin.immuno
                                   design = ~ condsp)
  dds.rp.deseq <- DESeq(dds.rp.deseq)
  res.rp.deseq <- results(dds.rp.deseq)
 
  resrp.deseq<-list(resr=res.rp.deseq)
  immuno.data.out.deseq.r[[i]] <- resrp.deseq #added [[i]] and referenced list in line prior
  
  #randomized with FP addition PD1
  dds.thp.deseq  <- DESeqDataSetFromMatrix(countData = datasp,
                                    colData = data.frame(condsp),
                                    design = ~ condsp)
  dds.thp.deseq <- DESeq(dds.thp.deseq)
  res.thp.deseq <- results(dds.thp.deseq)
  
  respp.deseq<-list(resp=res.thp.deseq)
  immuno.data.out.deseq.p[[i]] <- respp.deseq
 

  save(immuno.data.out.deseq.r, file="./analysis/immuno.data.r.deseq.Rda")
  save(immuno.data.out.deseq.p, file="./analysis/immuno.data.p.deseq.Rda")
  
}

#unpermuted PD1
dds.up.deseq  <- DESeqDataSetFromMatrix(countData = immuno.data,
                                 colData = immuno.conds,
                                 design = ~ conditions_p)
dds.up.deseq <- DESeq(dds.up.deseq)
res.up.deseq <- results(dds.up.deseq)

immuno.data.out.deseq.u<-list(resu=res.up.deseq)
save(immuno.data.out.deseq.u, file="./analysis/immuno.data.u.deseq.Rda") #}

