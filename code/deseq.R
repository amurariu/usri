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
data.out <- list() #check this, unsure of how this works

#brca
y_brca <- DGEList(counts=brca, group=factor(conditions_b))
keep_brca <- filterByExpr(y_brca)
y_brca <- y_brca[keep_brca,keep.lib.sizes=FALSE]
brca.data <- y_brca$counts #filtered base dataset
brca.data.out <- list() #check this, unsure of how this works

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
  
  
  data.iter.pd1.r <- list(desr=res.rp)
  pd1.data.out.r[[i]] <- data.iter.pd1.r
  save(pd1.data.out.r, file="./Documents/github/usri/analysis/pd1_data_out.r.Rda")
  
  data.iter.pd1.p <- list(desp=res.thp)
  pd1.data.out.p[[i]] <- data.iter.pd1.p
  save(pd1.data.out.p, file="./Documents/github/usri/analysis/pd1_data_out.p.Rda")
  
  data.iter.brca.r <- list(desr=res.rb)
  brca.data.out.r[[i]] <- data.iter.brca.r
  save(brca.data.out.r, file="./Documents/github/usri/analysis/brca_data_out.r.Rda")
  
  data.iter.brca.p <- list(desp=res.thb)
  brca.data.out.p[[i]] <- data.iter.brca.p
  save(brca.data.out.p, file="./Documents/github/usri/analysis/brca_data_out.p.Rda")
  
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

data.iter.pd1.u <- list(desp=res.up)
pd1.data.out.u[[i]] <- data.iter.pd1.u
save(pd1.data.out.u, file="./Documents/github/usri/analysis/pd1_data_out.u.Rda")

data.iter.brca.u <- list(desr=res.ub)
brca.data.out.u[[i]] <- data.iter.brca.u
save(brca.data.out.u, file="./Documents/github/usri/analysis/brca_data_out.u.Rda")


#start analysis from here
