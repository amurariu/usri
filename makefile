#makes the file .Rda if make all is called
#all: analysis/test2.Rda
DESeq: ../ext_analysis/immuno.data.deseq.Rda

edgeR: ../ext_analysis/immuno.data.edger.Rda

ALDEx2: ../ext_analysis/immuno.data.aldex2_5.out.Rda

ALDEx3: ../ext_analysis/immuno.data.aldex3.out.Rda

#ALDEx5: analysis/immuno.data.u.aldex5.Rda

#all: data_collection  

#rules to generate the deseq output files

../ext_analysis/immuno.data.deseq.Rda : code/deseq.R code/des.fun.R
	Rscript 'code/deseq.R' 'code/des.fun.R'

../ext_analysis/immuno.data.edger.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'

../ext_analysis/immuno.data.aldex2_5.out.Rda : code/aldex2_0.R code/ald2.fun.R
	Rscript 'code/aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/immuno.data.aldex3.out.Rda : code/aldex3_0.R code/ald3.fun.R
	Rscript 'code/aldex3_0.R' 'code/ald3.fun.R'

#analysis/immuno.data.u.aldex5.Rda : code/aldex_0.5.R
	#Rscript 'code/aldex_0.R'
	
#analysis/test2.Rda: code/brca_draft.R
#	Rscript 'code/brca_draft.R'

clean_DESeq:
	rm ../ext_analysis/immuno.data.deseq.Rda

clean_edgeR:
	rm ../ext_analysis/immuno.data.edger.Rda

clean_ALDEx2:
	#rm ../ext_analysis/immuno.data.aldex2_0.out.Rda
	#rm ../ext_analysis/immuno.data.aldex2_2.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_5.out.Rda
	

#	rm analysis/immuno.data.u.aldex2.Rda

#	rm analysis/immuno.data.u.aldex5.Rda

	
	
