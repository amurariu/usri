##### sanity checks

load("/Users/andreeamurariu/Documents/github/ext_analysis/immuno.data.deseq.Rda")
load("/Users/andreeamurariu/Documents/github/ext_analysis/immuno.data.edger.Rda")
load("/Users/andreeamurariu/Documents/github/ext_analysis/immuno.data.aldex2_0.out.Rda")
load("/Users/andreeamurariu/Documents/github/ext_analysis/immuno.data.aldex2_2.out.Rda")
load("/Users/andreeamurariu/Documents/github/ext_analysis/immuno.data.aldex2_5.out.Rda")

# "immuno.data_0.aldex2" "immuno.data_2.aldex2"
# "immuno.data_5.aldex2" "immuno.data.DESeq"    
# "immuno.data.edgeR"  

##### sanity check 1
## AM check for each vs edgeR
## only need to do for one random instance
# is the thin data identical in each analysis
plot(immuno.data.edgeR$thin.data[[11]]$coefmat, immuno.data.DESeq$thin.data[[11]]$coefmat)
plot(immuno.data.edgeR$thin.data[[11]]$coefmat, immuno.data_0.aldex2$thin.data[[11]]$coefmat)
plot(immuno.data.edgeR$thin.data[[11]]$coefmat, immuno.data_2.aldex2$thin.data[[11]]$coefmat)
plot(immuno.data.edgeR$thin.data[[11]]$coefmat, immuno.data_5.aldex2$thin.data[[11]]$coefmat) #### AM: does not look right, looks like cross instead of diagonal?

# do a correlation for each between samples - should be 1 
# AM: looks good for the first 3 above (not aldex 0.5)

# subsetting the 100 entries
# variable$dataset[[replicate]]$name
plot(immuno.data.edgeR$thin.data[[11]]$coefmat,immuno.data_5.aldex2$thin.data[[11]]$coefmat)

##### sanity check 2
# should be different within each group
# major density in the middle in a cross shape
# AM check a couple more in each analysis
plot(immuno.data.edgeR$thin.data[[23]]$coefmat,immuno.data.edgeR$thin.data[[54]]$coefmat, pch=19, cex=0.2, col=rgb(0,0,0,0.5))
plot(immuno.data.DESeq$thin.data[[26]]$coefmat,immuno.data.DESeq$thin.data[[55]]$coefmat, pch=19, cex=0.2, col=rgb(0,0,0,0.5))
plot(immuno.data_0.aldex2$thin.data[[80]]$coefmat,immuno.data_0.aldex2$thin.data[[73]]$coefmat, pch=19, cex=0.2, col=rgb(0,0,0,0.5))
plot(immuno.data_2.aldex2$thin.data[[40]]$coefmat,immuno.data_2.aldex2$thin.data[[20]]$coefmat, pch=19, cex=0.2, col=rgb(0,0,0,0.5))
plot(immuno.data_5.aldex2$thin.data[[15]]$coefmat,immuno.data_5.aldex2$thin.data[[44]]$coefmat, pch=19, cex=0.2, col=rgb(0,0,0,0.5))


###### sanity check 3
# similar metrics in each analysis should be similar but not identical
# logFC (E), log2FoldChange (D), diff.btw (A)
# logCPM (E), log(baseMean) (D), rab.all (A)
# FDR (E), padj (D), we.eBH (A), wi.eBH (A)
# AM check a few more

#for logFC, log2FoldChange, diff.btw - for r.data
#AM: seems fine, points centered around 0 for both axes
plot(immuno.data.edgeR$r.data[[11]]$logFC,immuno.data.DESeq$r.data[[11]]$log2FoldChange, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log="y")
plot(immuno.data.edgeR$r.data[[11]]$logFC,immuno.data_0.aldex2$r.data[[11]]$diff.btw, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log="y")
plot(immuno.data.edgeR$r.data[[11]]$logFC,immuno.data_2.aldex2$r.data[[11]]$diff.btw, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log="y")
plot(immuno.data.edgeR$r.data[[11]]$logFC,immuno.data_5.aldex2$r.data[[11]]$diff.btw, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log="y")

#for log CPM, log(BaseMean), rab.all 
#AM: also seems fine, get diagonal line for all 4 variations
plot(immuno.data.edgeR$r.data[[11]]$logCPM,immuno.data.DESeq$r.data[[11]]$baseMean, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log="y")
plot(immuno.data.edgeR$r.data[[11]]$logCPM,immuno.data_0.aldex2$r.data[[11]]$rab.all, pch=19, cex=0.2, col=rgb(0,0,0,0.5))
plot(immuno.data.edgeR$r.data[[11]]$logCPM,immuno.data_2.aldex2$r.data[[11]]$rab.all, pch=19, cex=0.2, col=rgb(0,0,0,0.5))
plot(immuno.data.edgeR$r.data[[11]]$logCPM,immuno.data_5.aldex2$r.data[[11]]$rab.all, pch=19, cex=0.2, col=rgb(0,0,0,0.5))

# only random data
# edgeR shows positives, ALDEx2 not
plot(immuno.data.edgeR$r.data[[30]]$FDR,immuno.data_0.aldex2$r.data[[30]]$we.eBH, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log='xy') #get odd graphs here for 30,32,33
abline(v=0.05)
abline(h=0.05)

plot(immuno.data.edgeR$r.data[[11]]$FDR,immuno.data_0.aldex2$r.data[[11]]$wi.eBH, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log='xy') #get odd graphs here for 30,32,33
abline(v=0.05)
abline(h=0.05)

# both false positives
plot(immuno.data.edgeR$r.data[[11]]$FDR,immuno.data.DESeq$r.data[[11]]$padj, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log='xy')
abline(v=0.05)
abline(h=0.05)


# random data + TP
# edgeR and ALDEx2 find similar TP
plot(immuno.data.edgeR$p.data[[11]]$FDR,immuno.data_5.aldex2$p.data[[11]]$we.eBH, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log='xy')
abline(v=0.05)
abline(h=0.05)

# both positives DESeq lower values than edgeR
plot(immuno.data.edgeR$p.data[[11]]$FDR,immuno.data.DESeq$p.data[[11]]$padj, pch=19, cex=0.2, col=rgb(0,0,0,0.5), log='xy')
abline(v=0.05)
abline(h=0.05)






