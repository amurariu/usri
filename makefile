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

ALDEx2_0_immuno: ../ext_analysis/immuno.data.aldex2_0.Rda
ALDEx2_0_brca: ../ext_analysis/brca.data.aldex2_0.Rda
ALDEx2_0_kirc: ../ext_analysis/kirc.data.aldex2_0.Rda
ALDEx2_0_lihc: ../ext_analysis/lihc.data.aldex2_0.Rda
ALDEx2_0_luad: ../ext_analysis/luad.data.aldex2_0.Rda
ALDEx2_0_prad: ../ext_analysis/prad.data.aldex2_0.Rda
ALDEx2_0_thca: ../ext_analysis/thca.data.aldex2_0.Rda

ALDEx2_1_immuno: ../ext_analysis/immuno.data.aldex2_1.Rda
ALDEx2_1_brca: ../ext_analysis/brca.data.aldex2_1.Rda
ALDEx2_1_kirc: ../ext_analysis/kirc.data.aldex2_1.Rda
ALDEx2_1_lihc: ../ext_analysis/lihc.data.aldex2_1.Rda
ALDEx2_1_luad: ../ext_analysis/luad.data.aldex2_1.Rda
ALDEx2_1_prad: ../ext_analysis/prad.data.aldex2_1.Rda
ALDEx2_1_thca: ../ext_analysis/thca.data.aldex2_1.Rda

ALDEx2_2_immuno: ../ext_analysis/immuno.data.aldex2_2.Rda
ALDEx2_2_brca: ../ext_analysis/brca.data.aldex2_2.Rda
ALDEx2_2_kirc: ../ext_analysis/kirc.data.aldex2_2.Rda
ALDEx2_2_lihc: ../ext_analysis/lihc.data.aldex2_2.Rda
ALDEx2_2_luad: ../ext_analysis/luad.data.aldex2_2.Rda
ALDEx2_2_prad: ../ext_analysis/prad.data.aldex2_2.Rda
ALDEx2_2_thca: ../ext_analysis/thca.data.aldex2_2.Rda

ALDEx2_5_immuno: ../ext_analysis/immuno.data.aldex2_5.Rda
ALDEx2_5_brca: ../ext_analysis/brca.data.aldex2_5.Rda
ALDEx2_5_kirc: ../ext_analysis/kirc.data.aldex2_5.Rda
ALDEx2_5_lihc: ../ext_analysis/lihc.data.aldex2_5.Rda
ALDEx2_5_luad: ../ext_analysis/luad.data.aldex2_5.Rda
ALDEx2_5_prad: ../ext_analysis/prad.data.aldex2_5.Rda
ALDEx2_5_thca: ../ext_analysis/thca.data.aldex2_5.Rda

limma_immuno: ../ext_analysis/immuno.data.limma.Rda
limma_brca: ../ext_analysis/brca.data.limma.Rda
limma_kirc: ../ext_analysis/kirc.data.limma.Rda
limma_lihc: ../ext_analysis/lihc.data.limma.Rda
limma_luad: ../ext_analysis/luad.data.limma.Rda
limma_prad: ../ext_analysis/prad.data.limma.Rda
limma_thca: ../ext_analysis/thca.data.limma.Rda

limma_immuno_asym1: ../ext_analysis/immuno.data.limma.asym1.Rda
limma_immuno_asym2: ../ext_analysis/immuno.data.limma.asym2.Rda
limma_immuno_asym3: ../ext_analysis/immuno.data.limma.asym3.Rda

limma_brca_asym1: ../ext_analysis/brca.data.limma.asym1.Rda
limma_brca_asym2: ../ext_analysis/brca.data.limma.asym2.Rda
limma_brca_asym3: ../ext_analysis/brca.data.limma.asym3.Rda

limma_kirc_asym1: ../ext_analysis/kirc.data.limma.asym1.Rda
limma_kirc_asym2: ../ext_analysis/kirc.data.limma.asym2.Rda
limma_kirc_asym3: ../ext_analysis/kirc.data.limma.asym3.Rda

limma_lihc_asym1: ../ext_analysis/lihc.data.limma.asym1.Rda
limma_lihc_asym2: ../ext_analysis/lihc.data.limma.asym2.Rda
limma_lihc_asym3: ../ext_analysis/lihc.data.limma.asym3.Rda

limma_luad_asym1: ../ext_analysis/luad.data.limma.asym1.Rda
limma_luad_asym2: ../ext_analysis/luad.data.limma.asym2.Rda
limma_luad_asym3: ../ext_analysis/luad.data.limma.asym3.Rda

limma_prad_asym1: ../ext_analysis/prad.data.limma.asym1.Rda
limma_prad_asym2: ../ext_analysis/prad.data.limma.asym2.Rda
limma_prad_asym3: ../ext_analysis/prad.data.limma.asym3.Rda

limma_thca_asym1: ../ext_analysis/thca.data.limma.asym1.Rda
limma_thca_asym2: ../ext_analysis/thca.data.limma.asym2.Rda
limma_thca_asym3: ../ext_analysis/thca.data.limma.asym3.Rda

#rules to generate the output files
### DESeq
../ext_analysis/immuno.data.deseq.Rda : code/immuno_deseq.R code/des.fun.R
	Rscript 'code/immuno_deseq.R' 'code/des.fun.R'

../ext_analysis/brca.data.deseq.Rda : code/brca_deseq.R code/des.fun.R
	Rscript 'code/brca_deseq.R' 'code/des.fun.R'
	
../ext_analysis/kirc.data.deseq.Rda : code/kirc_deseq.R code/des.fun.R
	Rscript 'code/kirc_deseq.R' 'code/des.fun.R'

../ext_analysis/lihc.data.deseq.Rda : code/lihc_deseq.R code/des.fun.R
	Rscript 'code/lihc_deseq.R' 'code/des.fun.R'

../ext_analysis/luad.data.deseq.Rda : code/luad_deseq.R code/des.fun.R
	Rscript 'code/luad_deseq.R' 'code/des.fun.R'

../ext_analysis/prad.data.deseq.Rda : code/prad_deseq.R code/des.fun.R
	Rscript 'code/prad_deseq.R' 'code/des.fun.R'

../ext_analysis/thca.data.deseq.Rda : code/thca_deseq.R code/des.fun.R
	Rscript 'code/thca_deseq.R' 'code/des.fun.R'


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


### aldex2_0
../ext_analysis/immuno.data.aldex2_0.Rda : code/immuno_aldex2_0.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/brca.data.aldex2_0.Rda : code/brca_aldex2_0.R code/ald2.fun.R
	Rscript 'code/brca_aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/kirc.data.aldex2_0.Rda : code/kirc_aldex2_0.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/lihc.data.aldex2_0.Rda : code/lihc_aldex2_0.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/luad.data.aldex2_0.Rda : code/luad_aldex2_0.R 
	Rscript 'code/luad_aldex2_0.R' 
	
../ext_analysis/prad.data.aldex2_0.Rda : code/prad_aldex2_0.R code/ald2.fun.R
	Rscript 'code/prad_aldex2_0.R' 'code/ald2.fun.R'
	
../ext_analysis/thca.data.aldex2_0.Rda : code/thca_aldex2_0.R code/ald2.fun.R
	Rscript 'code/thca_aldex2_0.R' 'code/ald2.fun.R'
	
### aldex2_1
../ext_analysis/immuno.data.aldex2_1.Rda : code/immuno_aldex2_1.R
	Rscript 'code/immuno_aldex2_1.R' 
	
../ext_analysis/brca.data.aldex2_1.Rda : code/brca_aldex2_1.R 
	Rscript 'code/brca_aldex2_1.R' 
	
../ext_analysis/kirc.data.aldex2_1.Rda : code/kirc_aldex2_1.R 
	Rscript 'code/kirc_aldex2_1.R'
	
../ext_analysis/lihc.data.aldex2_1.Rda : code/lihc_aldex2_1.R 
	Rscript 'code/lihc_aldex2_1.R' 
	
../ext_analysis/luad.data.aldex2_1.Rda : code/luad_aldex2_1.R 
	Rscript 'code/luad_aldex2_1.R' 
	
../ext_analysis/prad.data.aldex2_1.Rda : code/prad_aldex2_1.R 
	Rscript 'code/prad_aldex2_1.R' 
	
../ext_analysis/thca.data.aldex2_1.Rda : code/thca_aldex2_1.R 
	Rscript 'code/thca_aldex2_1.R' 
	
### aldex2_2
../ext_analysis/immuno.data.aldex2_2.Rda : code/immuno_aldex2_2.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2_2.R' 'code/ald2.fun.R'
	
../ext_analysis/brca.data.aldex2_2.Rda : code/brca_aldex2_2.R code/ald2.fun.R
	Rscript 'code/brca_aldex2_2.R' 'code/ald2.fun.R'
	
../ext_analysis/kirc.data.aldex2_2.Rda : code/kirc_aldex2_2.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2_2.R' 'code/ald2.fun.R'
	
../ext_analysis/lihc.data.aldex2_2.Rda : code/lihc_aldex2_2.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2_2.R' 'code/ald2.fun.R'
	
../ext_analysis/luad.data.aldex2_2.Rda : code/luad_aldex2_2.R code/ald2.fun.R
	Rscript 'code/luad_aldex2_2.R' 'code/ald2.fun.R'
	
../ext_analysis/prad.data.aldex2_2.Rda : code/prad_aldex2_2.R code/ald2.fun.R
	Rscript 'code/prad_aldex2_2.R' 'code/ald2.fun.R'
	
../ext_analysis/thca.data.aldex2_2.Rda : code/thca_aldex2_2.R code/ald2.fun.R
	Rscript 'code/thca_aldex2_2.R' 'code/ald2.fun.R'
	
### aldex2_5
../ext_analysis/immuno.data.aldex2_5.Rda : code/immuno_aldex2_5.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2_5.R' 'code/ald2.fun.R'
	
../ext_analysis/brca.data.aldex2_5.Rda : code/brca_aldex2_5.R code/ald2.fun.R
	Rscript 'code/brca_aldex2_5.R' 'code/ald2.fun.R'
	
../ext_analysis/kirc.data.aldex2_5.Rda : code/kirc_aldex2_5.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2_5.R' 'code/ald2.fun.R'
	
../ext_analysis/lihc.data.aldex2_5.Rda : code/lihc_aldex2_5.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2_5.R' 'code/ald2.fun.R'
	
../ext_analysis/luad.data.aldex2_5.Rda : code/luad_aldex2_5.R code/ald2.fun.R
	Rscript 'code/luad_aldex2_5.R' 'code/ald2.fun.R'
	
../ext_analysis/prad.data.aldex2_5.Rda : code/prad_aldex2_5.R code/ald2.fun.R
	Rscript 'code/prad_aldex2_5.R' 'code/ald2.fun.R'
	
../ext_analysis/thca.data.aldex2_5.Rda : code/thca_aldex2_5.R code/ald2.fun.R
	Rscript 'code/thca_aldex2_5.R' 'code/ald2.fun.R'
	

###limma 	
../ext_analysis/immuno.data.limma.Rda : code/immuno_limma.R code/lim.fun.R
	Rscript 'code/immuno_limma.R' 'code/lim.fun.R'
	
../ext_analysis/brca.data.limma.Rda : code/brca_limma.R code/lim.fun.R
	Rscript 'code/brca_limma.R' 'code/lim.fun.R'
	
../ext_analysis/kirc.data.limma.Rda : code/kirc_limma.R code/lim.fun.R
	Rscript 'code/kirc_limma.R' 'code/lim.fun.R'
	
../ext_analysis/lihc.data.limma.Rda : code/lihc_limma.R code/lim.fun.R
	Rscript 'code/lihc_limma.R' 'code/lim.fun.R'

../ext_analysis/luad.data.limma.Rda : code/luad_limma.R code/lim.fun.R
	Rscript 'code/luad_limma.R' 'code/lim.fun.R'
	
../ext_analysis/prad.data.limma.Rda : code/prad_limma.R code/lim.fun.R
	Rscript 'code/prad_limma.R' 'code/lim.fun.R'
	
../ext_analysis/thca.data.limma.Rda : code/thca_limma.R code/lim.fun.R
	Rscript 'code/thca_limma.R' 'code/lim.fun.R'

#asymmetrical limma files
../ext_analysis/immuno.data.limma.asym1.Rda : code/immuno_limma_asym1.R
	Rscript 'code/immuno_limma_asym1.R'
	
../ext_analysis/immuno.data.limma.asym2.Rda : code/immuno_limma_asym2.R
	Rscript 'code/immuno_limma_asym2.R'
	
../ext_analysis/immuno.data.limma.asym3.Rda : code/immuno_limma_asym3.R
	Rscript 'code/immuno_limma_asym3.R'
	

../ext_analysis/brca.data.limma.asym1.Rda : code/brca_limma_asym1.R
	Rscript 'code/brca_limma_asym1.R'
	
../ext_analysis/brca.data.limma.asym2.Rda : code/brca_limma_asym2.R
	Rscript 'code/brca_limma_asym2.R'
	
../ext_analysis/brca.data.limma.asym3.Rda : code/brca_limma_asym3.R
	Rscript 'code/brca_limma_asym3.R'

../ext_analysis/kirc.data.limma.asym1.Rda : code/kirc_limma_asym1.R
	Rscript 'code/kirc_limma_asym1.R'
	
../ext_analysis/kirc.data.limma.asym2.Rda : code/kirc_limma_asym2.R
	Rscript 'code/kirc_limma_asym2.R'
	
../ext_analysis/kirc.data.limma.asym3.Rda : code/kirc_limma_asym3.R
	Rscript 'code/kirc_limma_asym3.R'

../ext_analysis/lihc.data.limma.asym1.Rda : code/lihc_limma_asym1.R
	Rscript 'code/lihc_limma_asym1.R'
	
../ext_analysis/lihc.data.limma.asym2.Rda : code/lihc_limma_asym2.R
	Rscript 'code/lihc_limma_asym2.R'
	
../ext_analysis/lihc.data.limma.asym3.Rda : code/lihc_limma_asym3.R
	Rscript 'code/lihc_limma_asym3.R'

../ext_analysis/luad.data.limma.asym1.Rda : code/luad_limma_asym1.R
	Rscript 'code/luad_limma_asym1.R'
	
../ext_analysis/luad.data.limma.asym2.Rda : code/luad_limma_asym2.R
	Rscript 'code/luad_limma_asym2.R'
	
../ext_analysis/luad.data.limma.asym3.Rda : code/luad_limma_asym3.R
	Rscript 'code/luad_limma_asym3.R'

../ext_analysis/prad.data.limma.asym1.Rda : code/prad_limma_asym1.R
	Rscript 'code/prad_limma_asym1.R'
	
../ext_analysis/prad.data.limma.asym2.Rda : code/prad_limma_asym2.R
	Rscript 'code/prad_limma_asym2.R'
	
../ext_analysis/prad.data.limma.asym3.Rda : code/prad_limma_asym3.R
	Rscript 'code/prad_limma_asym3.R'

../ext_analysis/thca.data.limma.asym1.Rda : code/thca_limma_asym1.R
	Rscript 'code/thca_limma_asym1.R'
	
../ext_analysis/thca.data.limma.asym2.Rda : code/thca_limma_asym2.R
	Rscript 'code/thca_limma_asym2.R'
	
../ext_analysis/thca.data.limma.asym3.Rda : code/thca_limma_asym3.R
	Rscript 'code/thca_limma_asym3.R'


#DESeq clean files
clean_immuno_DESeq:
	rm ../ext_analysis/immuno.data.deseq.Rda
clean_brca_DESeq:
	rm ../ext_analysis/brca.data.deseq.Rda
clean_kirc_DESeq:
	rm ../ext_analysis/kirc.data.deseq.Rda
clean_lihc_DESeq:
	rm ../ext_analysis/lihc.data.deseq.Rda
clean_luad_DESeq:
	rm ../ext_analysis/luad.data.deseq.Rda
clean_prad_DESeq:
	rm ../ext_analysis/prad.data.deseq.Rda
clean_thca_DESeq:
	rm ../ext_analysis/thca.data.deseq.Rda

#edgeR clean files
clean_immuno_edgeR:
	rm ../ext_analysis/immuno.data.edger.Rda
clean_brca_edgeR:
	rm ../ext_analysis/brca.data.edger.Rda
clean_kirc_edgeR:
	rm ../ext_analysis/kirc.data.edger.Rda
clean_lihc_edgeR:
	rm ../ext_analysis/lihc.data.edger.Rda
clean_luad_edgeR:
	rm ../ext_analysis/luad.data.edger.Rda
clean_prad_edgeR:
	rm ../ext_analysis/prad.data.edger.Rda
clean_thca_edgeR:
	rm ../ext_analysis/thca.data.edger.Rda
	
#aldex clean files
clean_immuno_aldex2_0:
	rm ../ext_analysis/immuno.data.aldex2_0.Rda
clean_immuno_aldex2_1:
	rm ../ext_analysis/immuno.data.aldex2_1.Rda
clean_immuno_aldex2_2:
	rm ../ext_analysis/immuno.data.aldex2_2.Rda
clean_immuno_aldex2_5:
	rm ../ext_analysis/immuno.data.aldex2_5.Rda
clean_brca_aldex2_0:
	rm ../ext_analysis/brca.data.aldex2_0.Rda
clean_brca_aldex2_1:
	rm ../ext_analysis/brca.data.aldex2_1.Rda
clean_brca_aldex2_2:
	rm ../ext_analysis/brca.data.aldex2_2.Rda
clean_brca_aldex2_5:
	rm ../ext_analysis/brca.data.aldex2_5.Rda
clean_kirc_aldex2_0:
	rm ../ext_analysis/kirc.data.aldex2_0.Rda
clean_kirc_aldex2_1:
	rm ../ext_analysis/kirc.data.aldex2_1.Rda
clean_kirc_aldex2_2:
	rm ../ext_analysis/kirc.data.aldex2_2.Rda
clean_kirc_aldex2_5:
	rm ../ext_analysis/kirc.data.aldex2_5.Rda
clean_lihc_aldex2_0:
	rm ../ext_analysis/lihc.data.aldex2_0.Rda
clean_lihc_aldex2_1:
	rm ../ext_analysis/lihc.data.aldex2_1.Rda
clean_lihc_aldex2_2:
	rm ../ext_analysis/lihc.data.aldex2_2.Rda
clean_lihc_aldex2_5:
	rm ../ext_analysis/lihc.data.aldex2_5.Rda
clean_luad_aldex2_0:
	rm ../ext_analysis/luad.data.aldex2_0.Rda
clean_luad_aldex2_1:
	rm ../ext_analysis/luad.data.aldex2_1.Rda
clean_luad_aldex2_2:
	rm ../ext_analysis/luad.data.aldex2_2.Rda
clean_luad_aldex2_5:
	rm ../ext_analysis/luad.data.aldex2_5.Rda
clean_prad_aldex2_0:
	rm ../ext_analysis/prad.data.aldex2_0.Rda
clean_prad_aldex2_1:
	rm ../ext_analysis/prad.data.aldex2_1.Rda
clean_prad_aldex2_2:
	rm ../ext_analysis/prad.data.aldex2_2.Rda
clean_prad_aldex2_5:
	rm ../ext_analysis/prad.data.aldex2_5.Rda
clean_thca_aldex2_0:
	rm ../ext_analysis/thca.data.aldex2_0.Rda
clean_thca_aldex2_1:
	rm ../ext_analysis/thca.data.aldex2_1.Rda
clean_thca_aldex2_2:
	rm ../ext_analysis/thca.data.aldex2_2.Rda
clean_thca_aldex2_5:
	rm ../ext_analysis/thca.data.aldex2_5.Rda

#limma clean files
clean_immuno_limma:
	rm ../ext_analysis/immuno.data.limma.Rda
clean_brca_limma:
	rm ../ext_analysis/brca.data.limma.Rda
clean_kirc_limma:
	rm ../ext_analysis/kirc.data.limma.Rda
clean_lihc_limma:
	rm ../ext_analysis/lihc.data.limma.Rda
clean_luad_limma:
	rm ../ext_analysis/luad.data.limma.Rda
clean_prad_limma:
	rm ../ext_analysis/prad.data.limma.Rda
clean_thca_limma:
	rm ../ext_analysis/thca.data.limma.Rda

	
	
