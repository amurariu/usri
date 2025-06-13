#makes the file .Rda if make all is called
#all: analysis/test2.Rda
DESeq_immuno: ../ext_analysis/immuno.data.deseq.Rda
DESeq_brca: ../ext_analysis/brca.data.deseq.Rda

edgeR_immuno: ../ext_analysis/immuno.data.edger.Rda

ALDEx2_immuno: ../ext_analysis/immuno.data.aldex2_5.out.Rda

ALDEx3_immuno: ../ext_analysis/immuno.data.aldex3.out.Rda

limma_immuno: ../ext_analysis/immuno.data.limma.Rda


#rules to generate the output files

../ext_analysis/immuno.data.deseq.Rda : code/deseq.R code/des.fun.R
	Rscript 'code/deseq.R' 'code/des.fun.R'
	
../ext_analysis/brca.data.deseq.Rda : code/deseq.R code/des.fun.R
	Rscript 'code/deseq.R' 'code/des.fun.R'

../ext_analysis/immuno.data.edger.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'

../ext_analysis/immuno.data.aldex2_5.out.Rda : code/aldex2_0.R code/ald2.fun.R
	Rscript 'code/aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/immuno.data.aldex3.out.Rda : code/aldex3_0.R code/ald3.fun.R
	Rscript 'code/aldex3_0.R' 'code/ald3.fun.R'
	
../ext_analysis/immuno.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'


clean_DESeq:
	rm ../ext_analysis/immuno.data.deseq.Rda
	
clean_DESeq_brca:
	rm ../ext_analysis/brca.data.deseq.Rda


clean_edgeR:
	rm ../ext_analysis/immuno.data.edger.Rda

clean_ALDEx2:
	rm ../ext_analysis/immuno.data.aldex2_0.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_2.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_5.out.Rda
	
clean_limma:
	rm ../ext_analysis/immuno.data.limma.Rda
	


	
	
