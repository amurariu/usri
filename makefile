#makes the file .Rda if make all is called
#all: analysis/test2.Rda
DESeq: analysis/immuno.data.u.deseq.Rda

edgeR: analysis/immuno.data.u.edger.Rda

ALDEx: analysis/immuno.data.u.aldex0.Rda

ALDEx2: analysis/immuno.data.u.aldex2.Rda

ALDEx5: analysis/immuno.data.u.aldex5.Rda

#all: data_collection

#rules to generate the deseq output files

analysis/immuno.data.u.deseq.Rda : code/deseq.R
	Rscript 'code/deseq.R'

analysis/immuno.data.u.edger.Rda : code/edgeR.R
	Rscript 'code/edgeR.R'

analysis/immuno.data.u.aldex0.Rda : code/aldex_0.R
	Rscript 'code/aldex_0.R'
	
analysis/immuno.data.u.aldex2.Rda : code/aldex_0.2.R
	Rscript 'code/aldex_0.R'

analysis/immuno.data.u.aldex5.Rda : code/aldex_0.5.R
	Rscript 'code/aldex_0.R'
	
#analysis/test2.Rda: code/brca_draft.R
#	Rscript 'code/brca_draft.R'

clean all:
	rm analysis/immuno.data.u.deseq.Rda

#clean edgeR:
#	rm analysis/immuno.data.u.edger.Rda

#clean ALDEx:	
	rm analysis/immuno.data.u.aldex0.Rda

#clean ALDEx2:	
	rm analysis/immuno.data.u.aldex2.Rda

#clean ALDEx5:
	rm analysis/immuno.data.u.aldex5.Rda

	
	