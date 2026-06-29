Code to respond to reviewer 1 for PLoS Comp Bio

Two main issues: 16S not thinned so no standard of truth, and lack of single-cell bulk analysis.

16S analysis is in code/meta16S_aldex3.r

	This code solves the thinning problem. If it does not work, change line 26 to 500 features retained and re-try. This code retains those features that are closest to the median feature. 
	I have run this code with ALDEx3 and gamma at 0, and this example is in /data/meta16S_gamma0aldex3.Rda. I used the standard naming convention. 
	I also generated the confusion matrix in analysis/confusionMats//cm.meta16S.aldex3.0.Rda
	so, hopefully this can be used as a template to do the analysis with ALDEx2/3 to get the slope and intercept for 16S data
	
Pseudobulk analysis of single-cell data in code/pseudo.bulk_Seurat.R

	This code pulls single-cell pseudobulk data from the Seurat tutorial. This was a frustrating exercise because a) the Seurat manual is incomplete and requires a lot of reading between the lines b) I could not install the seurat-data R package using any described method. So instead, I cloned the github repository and load the package via devtools::load_all() as noted in the code. You should not need the R package as I have saved the needed data in data/ss.bulk.Rda. There is also one pre-computed ALDEx2 object.
	
	The data are the ifnb dataset in Seurat and are known single white blood cell types stimulated or not with some thing. The reference for the dataset is in the code. There are 103 control and 103 stimulated groups. 
	
	Line 92 and below shows the code needed to generate and extract the pseudobulk data.
	
Updated dispersion_abundance plot code in code/0_gg_fig_abundanceDispersion.R

	The layering needs fixing. The ss bulk data should be the first layer plotted.