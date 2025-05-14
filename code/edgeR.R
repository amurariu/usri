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
immuno.conds <- data.frame(conditions_p) 

#edgeR conditions for initial filtering
#PD1
y_pd1 <- DGEList(counts=immuno, group=factor(conditions_p))
keep_pd1 <- filterByExpr(y_pd1)
y_pd1 <- y_pd1[keep_pd1,keep.lib.sizes=FALSE]
immuno.data <- y_pd1$counts #filtered base dataset
immuno.data.out.edgeR.u <- list() 
immuno.data.out.edgeR.r <- list() 
immuno.data.out.edgeR.p <- list() 

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
  
  #edgeR analysis
  #PD1 setup
  group_p <- factor(condsp)
  design_p <- model.matrix(~group_p) #use data randomization from seqgendiff

  #randomized without FP addition PD1
  fit_rp <- glmQLFit(immuno.data,design_p) #uses original data (ie. no TP added)
  qlf_rp <- glmQLFTest(fit_rp,coef=2)
  edg.rp<-topTags(qlf_rp, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1)
 
  resrp.edgeR<-list(resu=edg.rp)
  immuno.data.out.edgeR.r[[i]] <- resrp.edgeR
  
  #randomized with FP addition PD1
  fit_pp <- glmQLFit(datasp,design_p)
  qlf_pp <- glmQLFTest(fit_pp,coef=2)
  edg.pp<-topTags(qlf_pp, n=nrow(datasp), adjust.method = "BH", sort.by = "none", p.value = 1)
  
  respp.edgeR<-list(resu=edg.pp)
  immuno.data.out.edgeR.p[[i]] <- respp.edgeR
  
}
  
#unpermuted PD1
group_up<-factor(conditions_p)
design_up <- model.matrix(~group_up)
fit_up <- glmQLFit(y_pd1,design_up)
qlf_up <- glmQLFTest(fit_up,coef=2)
edg.up<-topTags(qlf_up, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1) 

resup.edgeR<-list(resu=edg.up)
immuno.data.out.edgeR.u <- list(resup.edgeR)


#saving file
#PD1 save file
save(immuno.data.out.edgeR.u, file="./analysis/immuno.data.u.edger.Rda")
save(immuno.data.out.edgeR.r, file="./analysis/immuno.data.r.edger.Rda")
save(immuno.data.out.edgeR.p, file="./analysis/immuno.data.p.edger.Rda")
