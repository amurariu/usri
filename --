#makes the file .Rda if make all is called
#all: analysis/test2.Rda

nohub make &

#rules to generate the deseq output files
analysis/immuno.data.u.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'
	
#analysis/immuno.data.r.deseq.Rda: code/deseq.R
#	Rscript 'code/deseq.R'
	
#analysis/immuno.data.p.deseq.Rda: code/deseq.R
#	Rscript 'code/deseq.R'
	
#analysis/brca.data.u.deseq.Rda: code/deseq.R
#	Rscript 'code/deseq.R'

#analysis/brca.data.r.deseq.Rda: code/deseq.R
#	Rscript 'code/deseq.R'
	
#analysis/brca.data.p.deseq.Rda: code/deseq.R
#	Rscript 'code/deseq.R'
	
#analysis/test2.Rda: code/brca_draft.R
#	Rscript 'code/brca_draft.R'

clean:	
	#rm analysis/test2.Rda
	
	
	rm analysis/immuno.data.u.deseq.Rda
	
#	rm analysis/immuno.data.r.deseq.Rda
	
#	rm analysis/immuno.data.p.deseq.Rda
	
#	rm analysis/brca.data.u.deseq.Rda
	
#	rm analysis/brca.data.r.deseq.Rda
	
#	rm analysis/brca.data.p.deseq.Rda

	