#makes the file .Rda if make all is called
#all: analysis/test2.Rda
DESeq: analysis/immuno.data.Rda

edgeR: analysis/immuno.data.edger.out.Rda

ALDEx: analysis/immuno.data.aldex.out.Rda

#ALDEx2: analysis/immuno.data.u.aldex2.Rda

#ALDEx5: analysis/immuno.data.u.aldex5.Rda

#all: data_collection  

#rules to generate the deseq output files

analysis/immuno.data.Rda : code/deseq.R code/des.fun.R
	Rscript 'code/deseq.R' 'code/des.fun.R'

analysis/immuno.data.edger.out.Rda : code/edgeR.R code/edg.fun.R
	Rscript 'code/edgeR.R' 'code/edg.fun.R'

analysis/immuno.data.aldex.out.Rda : code/aldex_0.R code/ald.fun.R
	Rscript 'code/aldex_0.R' 'code/ald.fun.R'
	
#analysis/immuno.data.u.aldex2.Rda : code/aldex_0.2.R
	#Rscript 'code/aldex_0.R'

#analysis/immuno.data.u.aldex5.Rda : code/aldex_0.5.R
	#Rscript 'code/aldex_0.R'
	
#analysis/test2.Rda: code/brca_draft.R
#	Rscript 'code/brca_draft.R'

clean DESeq:
	rm analysis/immuno.data.Rda

clean edgeR:
	rm analysis/immuno.data.edger.out.Rda

clean ALDEx:
	rm analysis/immuno.data.aldex.out.Rda

#	rm analysis/immuno.data.u.aldex2.Rda

#	rm analysis/immuno.data.u.aldex5.Rda

	
	
