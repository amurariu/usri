#makes the file .Rda if make all is called
#all: analysis/test2.Rda

SCRIPT_DIR = code
ANALYSIS_DIR = analysis

SCRIPT = $(SCRIPT_DIR)/brca_draft.R
ANALYSIS = $(ANALYSIS_DIR)/test2.Rda


#target
all: $ANALYSIS

#generating output
$(OUTPUT): $SCRIPT
	mkdir -p $(ANALYSIS_DIR)
	Rscript $(SCRIPT) > $(ANALYSIS)

clean:
	rm -f $(OUTPUT)


#target is BRCA_Draft.Rda
#prereq is brca_draft.R
#recipe is to run the R script

#analysis/test2.Rda: 'code/brca_draft.R'
#	Rscript 'code/brca_draft.R'

#clean:	
#	rm analysis/test2.Rda
	