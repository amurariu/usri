#makes the file x.Rda if make all is called
all: BRCA_draft.Rda

#target is BRCA_Draft.Rda
#prereq is brca_draft.R
#recipe is to run the R script

BRCA_draft.Rda: brca_draft.R
	Rscript 'brca_draft.R'

#clean:	
#rm BRCA_Draft.Rda