library(ALDEx2, warn.conflicts=F)
library(seqgendiff, warn.conflicts=F)
library(edgeR, warn.conflicts=F)
library(DESeq2, warn.conflicts=F)

#immuno/PD1 dataset loading
raw_counts <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
meta <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm_metadata.txt'
immuno<-read.table(file=raw_counts, header = T, skip=35, sep='\t', row.names = 1)
m <- read.table(file=meta, header=F, row.names=1, sep='\t')
#establishing conditions for PD1
conditions_p <- rep("Pre", 109)
conditions_p[grep("_On",m)] <- "On"
immuno.conds <- data.frame(conditions_p) #changed conditions to conditions_p to be consistent across datasets

#edgeR conditions for initial filtering
#PD1
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset
immuno.data.out.aldex5.u <- list() 
immuno.data.out.aldex5.r <- list() 
immuno.data.out.aldex5.p <- list() 

#for loop
for (i in 1:100){
  #thin_2group adds rnorm noise to 5% of the transcripts, generates TPs in the dataset
  #generate thin_2group for each dataset as well as labelling for conditions and new dataset
  
  #PD1
  thin.immuno <- thin_2group(immuno.data, prop_null=0.95, alpha=0,
                             signal_fun = stats::rnorm, 
                             signal_params = list(mean = 0, sd = 2))
  condsp <- as.vector(thin.immuno$designmat)   # permuted and thinned conditions and data
  datasp <- thin.immuno$mat
  
  #ALDEX2 code added
  #randomized without FP addition PD1
  xrp.aldex5 <- aldex(immuno.data, conditions=condsp, gamma=0.5) #uses original dataset but permuted conditions
  
  #randomized with FP addition PD1
  xpp.aldex5 <- aldex(datasp, conditions=condsp, gamma=0.5) #uses new dataset with permuted conditions
  
}

#unpermuted datasets
#unpermuted PD1
xup.aldex5 <- aldex(immuno.data, conditions=immuno.conds, gamma=0.5)

#Save files here
#PD1 save file
immuno.data.out.aldex5.u<-list(resu=xup.aldex5)
save(immuno.data.out.aldex5.u, file="./analysis/immuno.data.u.aldex5.Rda")

immuno.data.out.aldex5.r<-list(resr=xrp.aldex5)
save(immuno.data.out.aldex5.r, file="./analysis/immuno.data.r.aldex5.Rda")

immuno.data.out.aldex5.p<-list(resp=xpp.aldex5)
save(immuno.data.out.aldex5.p, file="./analysis/immuno.data.p.aldex5.Rda")

