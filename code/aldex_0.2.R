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
immuno.data.out.aldex2.u <- list() 
immuno.data.out.aldex2.r <- list() 
immuno.data.out.aldex2.p <- list() 

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
  xrp.aldex2 <- aldex(immuno.data, conditions=condsp, gamma=0.2) #uses original dataset but permuted conditions
 
  resrp.aldex2<-list(resu=xrp.aldex2)
  immuno.data.out.r[[i]] <- list(resrp.aldex2)
  
  #randomized with FP addition PD1
  xpp.aldex2 <- aldex(datasp, conditions=condsp, gamma=0.2) #uses new dataset with permuted conditions
  
  respp.aldex2<-list(resu=xpp.aldex2)
  immuno.data.out.aldex2.p[[i]] <- list(respp.aldex2)
  
}

#unpermuted datasets
#unpermuted PD1
xup.aldex2 <- aldex(immuno.data, conditions=immuno.conds, gamma=0.2)
immuno.data.out.aldex2.u <- list(xup.aldex2)


#Save files here
#PD1 save file
save(immuno.data.out.aldex2.u, file="./analysis/immuno.data.u.aldex2.Rda")
save(immuno.data.out.aldex2.r, file="./analysis/immuno.data.r.aldex2.Rda")
save(immuno.data.out.aldex2.p, file="./analysis/immuno.data.p.aldex2.Rda")
