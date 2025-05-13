#makes the file .Rda if make all is called
#all: analysis/test2.Rda

#rules to generate the deseq output files
analysis/immuno.data.u.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'
	
analysis/immuno.data.r.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'
	
analysis/immuno.data.p.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'
	
analysis/brca.data.u.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'

analysis/brca.data.r.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'
	
analysis/brca.data.p.deseq.Rda: code/deseq.R
	Rscript 'code/deseq.R'
	
#rules to generate the edger output files
analysis/immuno.data.u.edger.Rda: code/edgeR.R
	Rscript 'code/edgeR.R'
	
analysis/immuno.data.r.edger.Rda: code/edgeR.R
	Rscript 'code/edgeR.R'
	
analysis/immuno.data.p.edger.Rda: code/edgeR.R
	Rscript 'code/edgeR.R'
	
analysis/brca.data.u.edger.Rda: code/edgeR.R
	Rscript 'code/edgeR.R'

analysis/brca.data.r.edger.Rda: code/edgeR.R
	Rscript 'code/edgeR.R'
	
analysis/brca.data.p.edger.Rda: code/edgeR.R
	Rscript 'code/edgeR.R'
	
	#rules to generate the aldex_0 output files
analysis/immuno.data.u.aldex0.Rda: code/aldex_0.R
	Rscript 'code/aldex_0.R'
	
analysis/immuno.data.r.aldex0.Rda: code/aldex_0.R
	Rscript 'code/aldex_0.R'
	
analysis/immuno.data.p.aldex0.Rda: code/aldex_0.R
	Rscript 'code/aldex_0.R'
	
analysis/brca.data.u.aldex0.Rda: code/aldex_0.R
	Rscript 'code/aldex_0.R'

analysis/brca.data.r.aldex0.Rda: code/aldex_0.R
	Rscript 'code/aldex_0.R'
	
analysis/brca.data.p.aldex0.Rda: code/aldex_0.R
	Rscript 'code/aldex_0.R'
	
#analysis/test2.Rda: code/brca_draft.R
#	Rscript 'code/brca_draft.R'

clean:	
	#rm analysis/test2.Rda
	
	
	rm analysis/immuno.data.u.deseq.Rda
	
	rm analysis/immuno.data.r.deseq.Rda
	
	rm analysis/immuno.data.p.deseq.Rda
	
	rm analysis/brca.data.u.deseq.Rda
	
	rm analysis/brca.data.r.deseq.Rda
	
	rm analysis/brca.data.p.deseq.Rda

	rm analysis/immuno.data.u.edger.Rda
	
	rm analysis/immuno.data.r.edger.Rda
	
	rm analysis/immuno.data.p.edger.Rda
	
	rm analysis/brca.data.u.edger.Rda
	
	rm analysis/brca.data.r.edger.Rda
	
	rm analysis/brca.data.p.edger.Rda
	
	rm analysis/immuno.data.u.aldex0.Rda
	
	rm analysis/immuno.data.r.aldex0.Rda
	
	rm analysis/immuno.data.p.aldex0.Rda
	
	rm analysis/brca.data.u.aldex0.Rda
	
	rm analysis/brca.data.r.aldex0.Rda
	
	rm analysis/brca.data.p.aldex0.Rda


	
	