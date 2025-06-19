#makes the file .Rda if make all is called
#all: analysis/test2.Rda
DESeq_immuno: ../ext_analysis/immuno.data.deseq.Rda
DESeq_brca: ../ext_analysis/brca.data.deseq.Rda
DESeq_kirc: ../ext_analysis/kirc.data.deseq.Rda
DESeq_lihc: ../ext_analysis/lihc.data.deseq.Rda
DESeq_luad: ../ext_analysis/luad.data.deseq.Rda
DESeq_prad: ../ext_analysis/prad.data.deseq.Rda
DESeq_thca: ../ext_analysis/thca.data.deseq.Rda

edgeR_immuno: ../ext_analysis/immuno.data.edger.Rda
edgeR_brca: ../ext_analysis/brca.data.edger.Rda
edgeR_kirc: ../ext_analysis/kirc.data.edger.Rda
edgeR_lihc: ../ext_analysis/lihc.data.edger.Rda
edgeR_luad: ../ext_analysis/luad.data.edger.Rda
edgeR_prad: ../ext_analysis/prad.data.edger.Rda
edgeR_thca: ../ext_analysis/thca.data.edger.Rda

ALDEx2_immuno: ../ext_analysis/immuno.data.aldex2_0.out.Rda
ALDEx2_brca: ../ext_analysis/brca.data.aldex2_0.out.Rda
ALDEx2_kirc: ../ext_analysis/kirc.data.aldex2_0.out.Rda
ALDEx2_lihc: ../ext_analysis/lihc.data.aldex2_0.out.Rda
ALDEx2_luad: ../ext_analysis/luad.data.aldex2_0.out.Rda
ALDEx2_prad: ../ext_analysis/prad.data.aldex2_0.out.Rda
ALDEx2_thca: ../ext_analysis/thca.data.aldex2_0.out.Rda

limma_immuno: ../ext_analysis/immuno.data.limma.Rda
limma_brca: ../ext_analysis/brca.data.limma.Rda
limma_kirc: ../ext_analysis/kirc.data.limma.Rda
limma_lihc: ../ext_analysis/lihc.data.limma.Rda
limma_luad: ../ext_analysis/luad.data.limma.Rda
limma_prad: ../ext_analysis/prad.data.limma.Rda
limma_thca: ../ext_analysis/thca.data.limma.Rda


#rules to generate the output files
### DESeq
../ext_analysis/immuno.data.deseq.Rda : code/immuno_deseq.R code/des.fun.R
	Rscript 'code/immuno_deseq.R' 'code/des.fun.R'

../ext_analysis/brca.data.deseq.Rda : code/immuno_deseq.R code/des.fun.R
	Rscript 'code/brca_deseq.R' 'code/des.fun.R'

### edgeR
../ext_analysis/immuno.data.edger.Rda : code/immuno_edgeR.R code/edg.fun.R
	Rscript 'code/immuno_edgeR.R' 'code/edg.fun.R'
	
../ext_analysis/brca.data.edger.Rda : code/brca_edgeR.R code/edg.fun.R
	Rscript 'code/brca_edgeR.R' 'code/edg.fun.R'
	
../ext_analysis/kirc.data.edger.Rda : code/kirc_edgeR.R code/edg.fun.R
	Rscript 'code/kirc_edgeR.R' 'code/edg.fun.R'

../ext_analysis/lihc.data.edger.Rda : code/lihc_edgeR.R code/edg.fun.R
	Rscript 'code/lihc_edgeR.R' 'code/edg.fun.R'
	
../ext_analysis/luad.data.edger.Rda : code/luad_edgeR.R code/edg.fun.R
	Rscript 'code/luad_edgeR.R' 'code/edg.fun.R'
	
../ext_analysis/prad.data.edger.Rda : code/prad_edgeR.R code/edg.fun.R
	Rscript 'code/prad_edgeR.R' 'code/edg.fun.R'

../ext_analysis/thca.data.edger.Rda : code/thca_edgeR.R code/edg.fun.R
	Rscript 'code/thca_edgeR.R' 'code/edg.fun.R'


### aldex2
../ext_analysis/immuno.data.aldex2_0.out.Rda : code/immuno_aldex2.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2.R' 'code/ald2.fun.R'
	
../ext_analysis/brca.data.aldex2_0.out.Rda : code/brca_aldex2.R code/ald2.fun.R
	Rscript 'code/brca_aldex2.R' 'code/ald2.fun.R'
	
	../ext_analysis/kirc.data.aldex2_0.out.Rda : code/kirc_aldex2.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2.R' 'code/ald2.fun.R'
	
../ext_analysis/lihc.data.aldex2_0.out.Rda : code/lihc_aldex2.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2.R' 'code/ald2.fun.R'
	
	../ext_analysis/luad.data.aldex2_0.out.Rda : code/luad_aldex2.R code/ald2.fun.R
	Rscript 'code/luad_aldex2.R' 'code/ald2.fun.R'
	
../ext_analysis/prad.data.aldex2_0.out.Rda : code/prad_aldex2.R code/ald2.fun.R
	Rscript 'code/prad_aldex2.R' 'code/ald2.fun.R'
	
../ext_analysis/thca.data.aldex2_0.out.Rda : code/thca_aldex2.R code/ald2.fun.R
	Rscript 'code/thca_aldex2.R' 'code/ald2.fun.R'
	


###limma 	
../ext_analysis/immuno.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'
	
../ext_analysis/brca.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'
	
../ext_analysis/kirc.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'
	
../ext_analysis/lihc.data.limma.Rda : code/limma.R code/lim.fun.R
	Rscript 'code/limma.R' 'code/lim.fun.R'


clean_DESeq:
	rm ../ext_analysis/immuno.data.deseq.Rda
	
clean_DESeq_brca:
	rm ../ext_analysis/brca.data.deseq.Rda

clean_edgeR:
	rm ../ext_analysis/immuno.data.edger.Rda
	
clean_edgeR_brca:
	rm ../ext_analysis/brca.data.edger.Rda

clean_ALDEx2:
	rm ../ext_analysis/immuno.data.aldex2_0.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_2.out.Rda
	rm ../ext_analysis/immuno.data.aldex2_5.out.Rda
	
clean_limma:
	rm ../ext_analysis/immuno.data.limma.Rda
	
clean_limma_brca:
	rm ../ext_analysis/brca.data.limma.Rda
	


	


	
	
