library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

#load datasets

#immuno/PD1 dataset loading
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p) #changed conditions to conditions_p to be consistent across datasets

#brca dataset loading
raw_counts<- "https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.normal-tumor.pair.rawCount.tsv"
con<-"https://raw.githubusercontent.com/amurariu/usri/main/data/TCGA-BRCA.conditions.tsv"

brca <- read.table(file=raw_counts, header=T, row.names=1, sep='\t')
conditions_b <- as.vector(unlist(read.table(file=con, sep='\t'))) #changed from brca.conds to conditions_b
brca.conds <- data.frame(conditions_b) #changed from conditions to brca.conds for consistency with PD1 dataset

#insert additional datasets + conditions here:


#edgeR conditions for initial filtering
#PD1
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset
imumuno.data.out.u <- list() 
imumuno.data.out.r <- list() 
imumuno.data.out.p <- list() 


#brca
y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset
brca.data.out.u <- list() 
brca.data.out.r <- list() 
brca.data.out.p <- list() 


#repeat adding edgeR conditions for each new dataset


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
  
  #BRCA
  thin.brca <- thin_2group(brca.data, prop_null=0.95, alpha=0,
                           signal_fun = stats::rnorm, signal_params = list(mean = 0, sd = 2))
  condsb <- as.vector(thin.brca$designmat)  # permuted and thinned conditions and data
  datasb <- thin.brca$mat #changed these from conds and datasp to condsb and datasb for consistency
  
 
  #DESeq2 analysis
  #randomized without FP addition PD1
  dds.rp  <- DESeqDataSetFromMatrix(countData = immuno.data,  #uses original data (no TP added)
                                   colData = data.frame(condsp), #uses data randomization order from thin.immuno
                                   design = ~ condsp)
  dds.rp <- DESeq(dds.rp)
  res.rp <- results(dds.rp)
  
  #randomized with FP addition PD1
  dds.thp  <- DESeqDataSetFromMatrix(countData = datasp,
                                    colData = data.frame(condsp),
                                    design = ~ condsp)
  dds.thp <- DESeq(dds.thp)
  res.thp <- results(dds.thp)
  
  
  #randomized without FP addition BRCA
  dds.rb  <- DESeqDataSetFromMatrix(countData = brca.data,  #uses original data (no TP added)
                                   colData = data.frame(condsb), #uses data randomization order from thin.brca
                                   design = ~ condsb)
  dds.rb <- DESeq(dds.rb)
  res.rb <- results(dds.rb)
  
  #randomized with FP addition BRCA
  dds.thb  <- DESeqDataSetFromMatrix(countData = thin.brca$mat, #data that includes TP
                                    colData = data.frame(condsb),
                                    design = ~ condsb)
  dds.thb <- DESeq(dds.thb)
  res.thb <- results(dds.thb)
  
}

#unpermuted datasets
#unpermuted PD1
dds.up  <- DESeqDataSetFromMatrix(countData = immuno.data,
                                 colData = immuno.conds,
                                 design = ~ conditions_p)
dds.up <- DESeq(dds.up)
res.up <- results(dds.up)

#unpermuted BRCA
dds.ub  <- DESeqDataSetFromMatrix(countData = brca.data,
                                 colData = brca.conds,
                                 design = ~ conditions_b)
dds.ub <- DESeq(dds.ub)
res.ub <- results(dds.ub)


#PD1 save file
desup<-list(desu=res.up)
immuno.data.out.u <- list(desup)
save(immuno.data.out.u, file="./Documents/github/usri/analysis/immuno.data.u.deseq.Rda")

desrp<-list(desr=res.rp)
immuno.data.out.r <- list(desrp)
save(immuno.data.out.r, file="./Documents/github/usri/analysis/immuno.data.r.deseq.Rda")

despp<-list(desp=res.thp)
immuno.data.out.p <- list(despnm vp)
save(immuno.data.out.p, file="./Documents/github/usri/analysis/immuno.data.u.deseq.Rda")

#BRCA save file
combined_b<-list(desu=res.ub)
brca.data.out <- list(combined_b)
save(brca.data.out, file="./Documents/github/usri/analysis/brca.data.deseq.Rda")

combined_b<-list(desu=res.ub, desr=res.rb, desp=res.thb)
brca.data.out <- list(combined_b)
save(brca.data.out, file="./Documents/github/usri/analysis/brca.data.deseq.Rda")

combined_b<-list(desu=res.ub, desr=res.rb, desp=res.thb)
brca.data.out <- list(combined_b)
save(brca.data.out, file="./Documents/github/usri/analysis/brca.data.deseq.Rda")

#start analysis from here
