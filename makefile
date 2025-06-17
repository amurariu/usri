#makes the file .Rda if make all is called
#all: analysis/test2.Rda
DESeq_immuno: ../ext_analysis/immuno.data.deseq.Rda
DESeq_brca: ../ext_analysis/brca.data.deseq.Rda
DESeq_kirc: ../ext_analysis/kirc.data.deseq.Rda
DESeq_lihc: ../ext_analysis/lihc.data.deseq.Rda

edgeR_immuno: ../ext_analysis/immuno.data.edger.Rda
edgeR_brca: ../ext_analysis/brca.data.edger.Rda
edgeR_kirc: ../ext_analysis/kirc.data.edger.Rda
edgeR_lihc: ../ext_analysis/lihc.data.edger.Rda

ALDEx2_immuno: ../ext_analysis/immuno.data.aldex2_0.out.Rda
ALDEx2_brca: ../ext_analysis/brca.data.aldex2_0.out.Rda
ALDEx2_kirc: ../ext_analysis/kirc.data.aldex2_0.out.Rda
ALDEx2_lihc: ../ext_analysis/lihc.data.aldex2_0.out.Rda


#ALDEx3_immuno: ../ext_analysis/immuno.data.aldex3.out.Rda

limma_immuno: ../ext_analysis/immuno.data.limma.Rda
limma_brca: ../ext_analysis/brca.data.limma.Rda
limma_kirc: ../ext_analysis/kirc.data.limma.Rda
limma_lihc: ../ext_analysis/lihc.data.limma.Rda


#rules to generate the output files
### DESeq
../ext_analysis/immuno.data.deseq.Rda : code/deseq.R code/des.fun.R
	Rscript 'code/deseq.R' 'code/des.fun.R'

../ext_analysis/brca.data.deseq.Rda : code/deseq.R code/des.fun.R
	Rscript 'code/deseq.R' 'code/des.fun.R'

### edgeR
../ext_analysis/immuno.data.edger.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'
	
../ext_analysis/brca.data.edger.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'
	
../ext_analysis/kirc.data.edger.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'

../ext_analysis/lihc.data.edger.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'


### aldex
../ext_analysis/immuno.data.aldex2_5.out.Rda : code/aldex2_0.R code/ald2.fun.R
	Rscript 'code/aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/brca.data.aldex2_0.out.Rda : code/aldex2_0.R code/ald2.fun.R
	Rscript 'code/aldex2_0.R' 'code/ald2.fun.R'
	
#../ext_analysis/immuno.data.aldex3.out.Rda : code/aldex3_0.R code/ald3.fun.R
	Rscript 'code/aldex3_0.R' 'code/ald3.fun.R'


###limma 	
../ext_analysis/immuno.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'
	
../ext_analysis/brca.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'
	
../ext_analysis/kirc.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'
	
../ext_analysis/lihc.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'


clean_DESeq:
	rm ../ext_analysis/immuno.data.deseq.Rda
	
clean_DESeq_brca:
	rm ../ext_analysis/brca.data.deseq.Rda

clean_edgeR:
	rm ../ext_analysis/immuno.data.edger.Rda
	
clean_edgeR_brca:
	rm ../ext_analysis/brca.data.edger.Rda

clean_ALDEx2:
	rm ../ext_analysis/immuno.data.aldex2_0.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_2.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_5.out.Rda
	
clean_limma:
	rm ../ext_analysis/immuno.data.limma.Rda
	
clean_limma_brca:
	rm ../ext_analysis/brca.data.limma.Rda
	


	


	
	
