#####
# USE THIS CODE AS EXAMPLE FOR REPLACEMENT FIG3
# GG
####

# read in immuno thinned immuno dataset aldex3 outputs 

rda <- paste0("https://github.com/amurariu/usri/raw/refs/heads/main/data/imm_gamma",
              0:5, "aldex", rep(c("2.Rda", "3.Rda"), each = 6))

# load directly from GitHub
for(i in rda){
  tf <- tempfile(fileext = ".Rda")
  download.file(url = i, destfile = tf, mode = "wb")
  load(tf)
  unlink(tf)
  rm(tf)
}

# get significant with gamma = 0 and 0.5
sig0 <- gamma0aldex3$p.val.adj < 0.01
sig5 <- gamma5aldex3$p.val.adj < 0.01

# get old-style dual cutoff significant
sig.old <- gamma0aldex3$p.val.adj < 0.01 & abs(gamma0aldex3$estimate) > 1
 
# dual cutoff plot
par(mfrow=c(2,1))
plot(gamma0aldex3$estimate, -log10(gamma0aldex3$p.val.adj), pch=19, col=rgb(0,0,0,0.2), cex=0.8,
  xlab='estimate', ylab='-log10(BH-p.value)', main='volcano')
points(gamma0aldex3$estimate[sig0], -log10(gamma0aldex3$p.val.adj[sig0]), pch=19, col=rgb(1,0,0,0.5), cex=0.8)
points(gamma0aldex3$estimate[sig.old], -log10(gamma0aldex3$p.val.adj[sig.old]), pch=19, col=rgb(0,0,1,0.8), cex=0.8)

abline(h=-log10(0.01), col='grey', lty=1, lwd=2)
abline(v=1, col='grey', lty=1, lwd=2)
abline(v=-1, col='grey', lty=1, lwd=2)

legend("topleft", legend=c(expression(paste("sig, not sig est > 1")), expression(paste("sig est > 1")), "not sig"), bg='white', pch=19, col=c("red", "blue", "grey"))

# This is an effect plot but with the X axis being estimate
plot(gamma0aldex3$estimate, gamma0aldex3$std.error, pch=19, col=rgb(0,0,0,0.2), cex=0.8, 
   xlab='estimate', ylab='SE estimate', main='effect')
points(gamma0aldex3$estimate[sig0], (gamma0aldex3$std.error[sig0]), pch=19, col=rgb(1,0,0,0.5), cex=0.8)
points(gamma0aldex3$estimate[sig.old], gamma0aldex3$std.error[sig.old], pch=19, col=rgb(0,0,1,0.8), cex=0.8)

abline(v=1, col='grey', lty=1, lwd=2)
abline(v=-1, col='grey', lty=1, lwd=2)

# This is the code for using gamma to control the FDR
# gamma cutoff
par(mfrow=c(2,1))
plot(gamma0aldex3$estimate, -log10(gamma0aldex3$p.val.adj), pch=19, col=rgb(0,0,0,0.2), cex=0.8,
  xlab='estimate', ylab='-log10(BH-p.value)', main='volcano')
points(gamma0aldex3$estimate[sig0], -log10(gamma0aldex3$p.val.adj[sig0]), pch=19, col=rgb(1,0,0,0.5), cex=0.8)
points(gamma0aldex3$estimate[sig5], -log10(gamma0aldex3$p.val.adj[sig5]), pch=19, col=rgb(0,0,1,0.8), cex=0.8)

abline(h=-log10(0.01), col='grey', lty=1, lwd=2)
#abline(v=1, col='grey', lty=1, lwd=2)
#abline(v=-1, col='grey', lty=1, lwd=2)

legend("topleft", legend=c(expression(paste("sig ", gamma, "=0, not sig ", gamma, "=0.5")), expression(paste("sig ", gamma, "=0.5")), "not sig"), bg='white', pch=19, col=c("red", "blue", "grey"))



# effect plot, again estimate on the X axis
plot(gamma0aldex3$estimate, gamma0aldex3$std.error, pch=19, col=rgb(0,0,0,0.2), cex=0.8, 
   xlab='estimate', ylab='SE estimate', main='effect')
points(gamma0aldex3$estimate[sig0], (gamma0aldex3$std.error[sig0]), pch=19, col=rgb(1,0,0,0.5), cex=0.8)
points(gamma0aldex3$estimate[sig5], gamma0aldex3$std.error[sig5], pch=19, col=rgb(0,0,1,0.8), cex=0.8)
# approximate isarithm for gamma=0
abline(0,-0.3, col='grey', lty=2, lwd=2)
abline(0,0.3, col='grey', lty=2, lwd=2)
# approximate isarithm for gamma=0.5
abline(-.3,-0.3,col='grey', lty=3, lwd=2)
abline(-.3,0.3,col='grey', lty=3, lwd=2)


##### ALDEx3 on vaginal dataset
### WARNING: paths in setup.R need changing along with hard-coded paths below
### REQUIRES ALDEx3 which can be obtained from CRAN
library(ALDEx3)

devtools::load_all('~/Documents/0_git/CoDaSeq/CoDaSeq')
path.to.github <- "~/Documents/0_git/projects/dossantos2024study/"

# set path in setup.R to gg
source(paste(path.to.github, "code/setup.R", sep = ""))

# load in vNumber -> KEGG pathway lookup table from VIRGO
path.table <- read.table(paste(locn,'1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""), 
                         sep="\t", header=T, row.names=1, fill=TRUE)

# load in vector containing cst info from london/europe species heatmap
load(paste(path.to.github, 'Rdata/hm.metadata.Rda', sep = ""))

# load in london/europe heatmap column colour bar list
load(paste(path.to.github, "Rdata/hm.column.cols.Rda",sep = ""))

# remove two samples from the filtered london/europe feature table aggregated by
# K0 number (classed as BV but almost no BV organisms):
#   -  v.001A: close to 100% L. gasseri with practically no BV organisms
#   -  v.019A: around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
ko.both<-ko.both[,-c(9,22)]

# remove non-bacterial KOs from the K number-aggregated, filtered feature table
# (these were discovered during curation of 'Unknown' pathways)
#    - K03364: Eukaryotic cell division cycle 20-like protein 1
#    - K13963: Serpin B (eukaryotic serine protease inhibitor)
#    - K01173: Mitochondrial endonuclease G
#    - K12373: Human lysosomal hexosaminidase
#    - K14327: Eukaryotic regulator of nonsense transcripts 2
#    - K00863: Human triose/dihydroxyacetone kinase
#    - K00599: Eukaryotic tRNA N(3)-methylcytidine methyltransferase
#    - K13993: Human HSP20
#    - K00811: Chloroplastic aspartate aminotransferase
#    - K03260: Eukaryotic translation initiation factor 4G
#    - K00985: Enterovirus RNA-directed RNA polymerase

ko.both <- ko.both[-which(grepl(paste("K03364","K13963","K01173","K12373",
                                      "K14327","K00863","K00599","K13993",
                                      "K00811","K03260","K00985", sep = "|"),
                                rownames(ko.both))),]
# clr.sm modified to accept a constant for the offset of the estimate
clr.const.sm <- function(X, logComp, const=0, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  logScale <- -colMeans(logComp, dims=1)

  tmp <- P*nsample
  LambdaScale <- matrix(rnorm(tmp,const,gamma), P, nsample)
  logScale <- logScale + t(X)%*% LambdaScale
  return(logScale)
}


ko.conds <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 

X <- data.frame(ko.conds)

# generate ALDEx3 outputs for plotting
vag.clr <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=clr.sm, gamma=0)
vag.clr5 <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=clr.const.sm, const=-3, gamma=0.5)
vag.tss <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=tss.sm, gamma=0)
vag.tss5 <- aldex(ko.both, ~ko.conds, X, nsample=128, scale=tss.sm, gamma=0.5)

# shows that tss is the same as clr with offset of -3
par(mfrow=c(2,2))
aldex.plot(vag.clr, contrast='ko.conds', main='clr g=0', plot='eff', cex=0.6)
#abline(h=3, lty=2, col='grey')
aldex.plot(vag.tss, contrast='ko.conds', main='tss g=0', plot='eff', cex=0.6)
aldex.plot(vag.clr5, contrast='ko.conds', main='clr c=-3, g=0.5', plot='eff', cex=0.6)
aldex.plot(vag.tss5, contrast='ko.conds', main='tss g=0.5', plot='eff', cex=0.6)

sum.clr5 <- summary(vag.clr5)
sum.tss5 <- summary(vag.tss5)


