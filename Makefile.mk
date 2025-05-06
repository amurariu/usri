#makes the file x.Rda if make all is called
all: BRCA_Draft.Rda

#target is BRCA_Draft.Rda
#prereq is brca_draft.R
#recipe is to run the R script

x.Rda: brca_draft.R
	Rscript 'brca_draft.R'

clean:
	rm BRCA_Draft.Rda