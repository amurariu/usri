#makes the file .Rda if make all is called
all: analysis/test2.Rda

#target is BRCA_Draft.Rda
#prereq is brca_draft.R
#recipe is to run the R script

analysis/test2.Rda: code/brca_draft.R
	Rscript 'code/brca_draft.R'

clean:	
	rm analysis/test2.Rda
	