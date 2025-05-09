library(ALDEx2, warn.conflicts=F)
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
  

  #ALDEX2 code added
  #randomized without FP addition PD1
  xrp <- aldex(immuno.data, conditions=condsp, gamma=1e-3) #uses original dataset but permuted conditions

  #randomized with FP addition PD1
  xpp <- aldex(datasp, conditions=condsp, gamma=1e-3) #uses new dataset with permuted conditions
  
  #randomized without FP addition BRCA
  xrb <- aldex(brca.data, conditions=condsb, gamma=1e-3) #uses original dataset but permuted conditions
  
  #randomized with FP addition PD1
  xpb <- aldex(datasb, conditions=condsb, gamma=1e-3) #uses new dataset with permuted conditions
  
  #add code to save each file separately
  
}

#add unpermuted datasets