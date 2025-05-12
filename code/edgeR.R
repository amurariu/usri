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
immuno.data.out <- list() #check this, unsure of how this works

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
  
  #edgeR analysis
  #PD1 setup
  group_p <- factor(condsp)
  design_p <- model.matrix(~group_p) #use data randomization from seqgendiff

  #randomized without FP addition PD1
  fit_rp <- glmQLFit(immuno.data,design_p) #uses original data (ie. no TP added)
  qlf_rp <- glmQLFTest(fit_rp,coef=2)
  edg.rp<-topTags(qlf_rp, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1)
 
  resrp<-list(resu=edg.rp)
  immuno.data.out.r[[i]] <- resrp
  
  #randomized with FP addition PD1
  fit_pp <- glmQLFit(datasp,design_p)
  qlf_pp <- glmQLFTest(fit_pp,coef=2)
  edg.pp<-topTags(qlf_pp, n=nrow(datasp), adjust.method = "BH", sort.by = "none", p.value = 1)
  
  respp<-list(resu=edg.pp)
  immuno.data.out.p[[i]] <- respp
  
  #BRCA setup
  group_b <- factor(condsb)
  design_b <- model.matrix(~group_b) #use data randomization from seqgendiff
  
  #randomized without FP addition BRCA
  fit_rb <- glmQLFit(brca.data,design_b) #uses original data (ie. no TP added)
  qlf_rb <- glmQLFTest(fit_rb,coef=2)
  edg.rb<-topTags(qlf_rb, n=nrow(brca.data), adjust.method = "BH", sort.by = "none", p.value = 1)
  
  resrb<-list(resu=edg.rb)
  brca.data.out.r[[i]] <- resrb
  
  #randomized with FP addition BRCA
  fit_pb <- glmQLFit(datasb,design_b) #data with TPs
  qlf_pb <- glmQLFTest(fit_pb,coef=2)
  edg.pb<-topTags(qlf_pb, n=nrow(datasb), adjust.method = "BH", sort.by = "none", p.value = 1)
 
  respb<-list(resu=edg.pb)
  brca.data.out.p[[i]] <- respb
  
  #add code to save these each as separate files 
}
  

#unpermuted PD1
group_up<-factor(conditions_p)
design_up <- model.matrix(~group_up)
fit_up <- glmQLFit(y_pd1,design_up)
qlf_up <- glmQLFTest(fit_up,coef=2)
edg.up<-topTags(qlf_up, n=nrow(immuno.data), adjust.method = "BH", sort.by = "none", p.value = 1) 

resup<-list(resu=edg.up)
immuno.data.out.u <- list(resup)

#unpermuted BRCA
group_ub<-factor(conditions_b)
design_ub <- model.matrix(~group_ub)
fit_ub <- glmQLFit(y_brca,design_ub)
qlf_ub <- glmQLFTest(fit_ub,coef=2)
edg.ub<-topTags(qlf_ub, n=nrow(brca.data), adjust.method = "BH", sort.by = "none", p.value = 1)

resub<-list(resu=edg.ub)
brca.data.out.u <- list(resub)

#saving file
#PD1 save file
save(immuno.data.out.u, file="./Documents/github/usri/analysis/immuno.data.u.edger.Rda")
save(immuno.data.out.r, file="./Documents/github/usri/analysis/immuno.data.r.edger.Rda")
save(immuno.data.out.p, file="./Documents/github/usri/analysis/immuno.data.p.edger.Rda")

#BRCA save file
save(brca.data.out.u, file="./Documents/github/usri/analysis/brca.data.u.edger.Rda")
save(brca.data.out.r, file="./Documents/github/usri/analysis/brca.data.r.edger.Rda")
save(brca.data.out.p, file="./Documents/github/usri/analysis/brca.data.p.edger.Rda")
