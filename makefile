#makes the file .Rda if make all is called
#all: analysis/test2.Rda

out_path="/Volumes/data2/andreea"
DESeq_immuno: $out_path/ext_analysis/immuno.data.deseq.Rda
DESeq_brca: $out_path/ext_analysis/brca.data.deseq.Rda
DESeq_kirc: $out_path/ext_analysis/kirc.data.deseq.Rda
DESeq_lihc: $out_path/ext_analysis/lihc.data.deseq.Rda
DESeq_luad: $out_path/ext_analysis/luad.data.deseq.Rda
DESeq_prad: $out_path/ext_analysis/prad.data.deseq.Rda
DESeq_thca: $out_path/ext_analysis/thca.data.deseq.Rda

edgeR_immuno: $out_path/ext_analysis/immuno.data.edger.Rda
edgeR_brca: $out_path/ext_analysis/brca.data.edger.Rda
edgeR_kirc: $out_path/ext_analysis/kirc.data.edger.Rda
edgeR_lihc: $out_path/ext_analysis/lihc.data.edger.Rda
edgeR_luad: $out_path/ext_analysis/luad.data.edger.Rda
edgeR_prad: $out_path/ext_analysis/prad.data.edger.Rda
edgeR_thca: $out_path/ext_analysis/thca.data.edger.Rda

ALDEx2_0_immuno: $out_path/ext_analysis/immuno.data.aldex2_0.Rda
ALDEx2_0_brca: $out_path/ext_analysis/brca.data.aldex2_0.Rda
ALDEx2_0_kirc: $out_path/ext_analysis/kirc.data.aldex2_0.Rda
ALDEx2_0_lihc: $out_path/ext_analysis/lihc.data.aldex2_0.Rda
ALDEx2_0_luad: $out_path/ext_analysis/luad.data.aldex2_0.Rda
ALDEx2_0_prad: $out_path/ext_analysis/prad.data.aldex2_0.Rda
ALDEx2_0_thca: $out_path/ext_analysis/thca.data.aldex2_0.Rda

ALDEx2_1_immuno: $out_path/ext_analysis/immuno.data.aldex2_1.Rda
ALDEx2_1_brca: $out_path/ext_analysis/brca.data.aldex2_1.Rda
ALDEx2_1_kirc: $out_path/ext_analysis/kirc.data.aldex2_1.Rda
ALDEx2_1_lihc: $out_path/ext_analysis/lihc.data.aldex2_1.Rda
ALDEx2_1_luad: $out_path/ext_analysis/luad.data.aldex2_1.Rda
ALDEx2_1_prad: $out_path/ext_analysis/prad.data.aldex2_1.Rda
ALDEx2_1_thca: $out_path/ext_analysis/thca.data.aldex2_1.Rda

ALDEx2_2_immuno: $out_path/ext_analysis/immuno.data.aldex2_2.Rda
ALDEx2_2_brca: $out_path/ext_analysis/brca.data.aldex2_2.Rda
ALDEx2_2_kirc: $out_path/ext_analysis/kirc.data.aldex2_2.Rda
ALDEx2_2_lihc: $out_path/ext_analysis/lihc.data.aldex2_2.Rda
ALDEx2_2_luad: $out_path/ext_analysis/luad.data.aldex2_2.Rda
ALDEx2_2_prad: $out_path/ext_analysis/prad.data.aldex2_2.Rda
ALDEx2_2_thca: $out_path/ext_analysis/thca.data.aldex2_2.Rda

ALDEx2_5_immuno: $out_path/ext_analysis/immuno.data.aldex2_5.Rda
ALDEx2_5_brca: $out_path/ext_analysis/brca.data.aldex2_5.Rda
ALDEx2_5_kirc: $out_path/ext_analysis/kirc.data.aldex2_5.Rda
ALDEx2_5_lihc: $out_path/ext_analysis/lihc.data.aldex2_5.Rda
ALDEx2_5_luad: $out_path/ext_analysis/luad.data.aldex2_5.Rda
ALDEx2_5_prad: $out_path/ext_analysis/prad.data.aldex2_5.Rda
ALDEx2_5_thca: $out_path/ext_analysis/thca.data.aldex2_5.Rda

limma_immuno: $out_path/ext_analysis/immuno.data.limma.Rda
limma_brca: $out_path/ext_analysis/brca.data.limma.Rda
limma_kirc: $out_path/ext_analysis/kirc.data.limma.Rda
limma_lihc: $out_path/ext_analysis/lihc.data.limma.Rda
limma_luad: $out_path/ext_analysis/luad.data.limma.Rda
limma_prad: $out_path/ext_analysis/prad.data.limma.Rda
limma_thca: $out_path/ext_analysis/thca.data.limma.Rda

limma_immuno_asym1: $out_path/ext_analysis/immuno.data.limma.asym1.Rda
limma_immuno_asym2: $out_path/ext_analysis/immuno.data.limma.asym2.Rda
limma_immuno_asym3: $out_path/ext_analysis/immuno.data.limma.asym3.Rda

limma_brca_asym1: $out_path/ext_analysis/brca.data.limma.asym1.Rda
limma_brca_asym2: $out_path/ext_analysis/brca.data.limma.asym2.Rda
limma_brca_asym3: $out_path/ext_analysis/brca.data.limma.asym3.Rda

limma_kirc_asym1: $out_path/ext_analysis/kirc.data.limma.asym1.Rda
limma_kirc_asym2: $out_path/ext_analysis/kirc.data.limma.asym2.Rda
limma_kirc_asym3: $out_path/ext_analysis/kirc.data.limma.asym3.Rda

limma_lihc_asym1: $out_path/ext_analysis/lihc.data.limma.asym1.Rda
limma_lihc_asym2: $out_path/ext_analysis/lihc.data.limma.asym2.Rda
limma_lihc_asym3: $out_path/ext_analysis/lihc.data.limma.asym3.Rda

limma_luad_asym1: $out_path/ext_analysis/luad.data.limma.asym1.Rda
limma_luad_asym2: $out_path/ext_analysis/luad.data.limma.asym2.Rda
limma_luad_asym3: $out_path/ext_analysis/luad.data.limma.asym3.Rda

limma_prad_asym1: $out_path/ext_analysis/prad.data.limma.asym1.Rda
limma_prad_asym2: $out_path/ext_analysis/prad.data.limma.asym2.Rda
limma_prad_asym3: $out_path/ext_analysis/prad.data.limma.asym3.Rda

limma_thca_asym1: $out_path/ext_analysis/thca.data.limma.asym1.Rda
limma_thca_asym2: $out_path/ext_analysis/thca.data.limma.asym2.Rda
limma_thca_asym3: $out_path/ext_analysis/thca.data.limma.asym3.Rda

ALDEx3_0_immuno: $out_path/ext_analysis/immuno.data.aldex3_0.Rda
ALDEx3_0_brca: $out_path/ext_analysis/brca.data.aldex3_0.Rda
ALDEx3_0_kirc: $out_path/ext_analysis/kirc.data.aldex3_0.Rda
ALDEx3_0_lihc: $out_path/ext_analysis/lihc.data.aldex3_0.Rda
ALDEx3_0_luad: $out_path/ext_analysis/luad.data.aldex3_0.Rda
ALDEx3_0_prad: $out_path/ext_analysis/prad.data.aldex3_0.Rda
ALDEx3_0_thca: $out_path/ext_analysis/thca.data.aldex3_0.Rda

ALDEx3_1_immuno: $out_path/ext_analysis/immuno.data.aldex3_1.Rda
ALDEx3_1_brca: $out_path/ext_analysis/brca.data.aldex3_1.Rda
ALDEx3_1_kirc: $out_path/ext_analysis/kirc.data.aldex3_1.Rda
ALDEx3_1_lihc: $out_path/ext_analysis/lihc.data.aldex3_1.Rda
ALDEx3_1_luad: $out_path/ext_analysis/luad.data.aldex3_1.Rda
ALDEx3_1_prad: $out_path/ext_analysis/prad.data.aldex3_1.Rda
ALDEx3_1_thca: $out_path/ext_analysis/thca.data.aldex3_1.Rda

ALDEx3_2_immuno: $out_path/ext_analysis/immuno.data.aldex3_2.Rda
ALDEx3_2_brca: $out_path/ext_analysis/brca.data.aldex3_2.Rda
ALDEx3_2_kirc: $out_path/ext_analysis/kirc.data.aldex3_2.Rda
ALDEx3_2_lihc: $out_path/ext_analysis/lihc.data.aldex3_2.Rda
ALDEx3_2_luad: $out_path/ext_analysis/luad.data.aldex3_2.Rda
ALDEx3_2_prad: $out_path/ext_analysis/prad.data.aldex3_2.Rda
ALDEx3_2_thca: $out_path/ext_analysis/thca.data.aldex3_2.Rda

ALDEx3_5_immuno: $out_path/ext_analysis/immuno.data.aldex3_5.Rda
ALDEx3_5_brca: $out_path/ext_analysis/brca.data.aldex3_5.Rda
ALDEx3_5_kirc: $out_path/ext_analysis/kirc.data.aldex3_5.Rda
ALDEx3_5_lihc: $out_path/ext_analysis/lihc.data.aldex3_5.Rda
ALDEx3_5_luad: $out_path/ext_analysis/luad.data.aldex3_5.Rda
ALDEx3_5_prad: $out_path/ext_analysis/prad.data.aldex3_5.Rda
ALDEx3_5_thca: $out_path/ext_analysis/thca.data.aldex3_5.Rda

limma_immuno_asym0_90: $out_path/ext_analysis/immuno.data.limma.asym0.prop90.Rda
limma_immuno_asym0_80: $out_path/ext_analysis/immuno.data.limma.asym0.prop80.Rda
limma_immuno_asym0_70: $out_path/ext_analysis/immuno.data.limma.asym0.prop70.Rda

limma_immuno_asym1_90: $out_path/ext_analysis/immuno.data.limma.asym1.prop90.Rda
limma_immuno_asym1_80: $out_path/ext_analysis/immuno.data.limma.asym1.prop80.Rda
limma_immuno_asym1_70: $out_path/ext_analysis/immuno.data.limma.asym1.prop70.Rda

limma_brca_asym0_90: $out_path/ext_analysis/brca.data.limma.asym0.prop90.Rda
limma_brca_asym0_80: $out_path/ext_analysis/brca.data.limma.asym0.prop80.Rda
limma_brca_asym0_70: $out_path/ext_analysis/brca.data.limma.asym0.prop70.Rda

limma_brca_asym1_90: $out_path/ext_analysis/brca.data.limma.asym1.prop90.Rda
limma_brca_asym1_80: $out_path/ext_analysis/brca.data.limma.asym1.prop80.Rda
limma_brca_asym1_70: $out_path/ext_analysis/brca.data.limma.asym1.prop70.Rda

limma_kirc_asym0_90: $out_path/ext_analysis/kirc.data.limma.asym0.prop90.Rda
limma_kirc_asym0_80: $out_path/ext_analysis/kirc.data.limma.asym0.prop80.Rda
limma_kirc_asym0_70: $out_path/ext_analysis/kirc.data.limma.asym0.prop70.Rda

limma_kirc_asym1_90: $out_path/ext_analysis/kirc.data.limma.asym1.prop90.Rda
limma_kirc_asym1_80: $out_path/ext_analysis/kirc.data.limma.asym1.prop80.Rda
limma_kirc_asym1_70: $out_path/ext_analysis/kirc.data.limma.asym1.prop70.Rda

limma_lihc_asym0_90: $out_path/ext_analysis/lihc.data.limma.asym0.prop90.Rda
limma_lihc_asym0_80: $out_path/ext_analysis/lihc.data.limma.asym0.prop80.Rda
limma_lihc_asym0_70: $out_path/ext_analysis/lihc.data.limma.asym0.prop70.Rda

limma_lihc_asym1_90: $out_path/ext_analysis/lihc.data.limma.asym1.prop90.Rda
limma_lihc_asym1_80: $out_path/ext_analysis/lihc.data.limma.asym1.prop80.Rda
limma_lihc_asym1_70: $out_path/ext_analysis/lihc.data.limma.asym1.prop70.Rda

limma_luad_asym0_90: $out_path/ext_analysis/luad.data.limma.asym0.prop90.Rda
limma_luad_asym0_80: $out_path/ext_analysis/luad.data.limma.asym0.prop80.Rda
limma_luad_asym0_70: $out_path/ext_analysis/luad.data.limma.asym0.prop70.Rda

limma_luad_asym1_90: $out_path/ext_analysis/luad.data.limma.asym1.prop90.Rda
limma_luad_asym1_80: $out_path/ext_analysis/luad.data.limma.asym1.prop80.Rda
limma_luad_asym1_70: $out_path/ext_analysis/luad.data.limma.asym1.prop70.Rda

limma_prad_asym0_90: $out_path/ext_analysis/prad.data.limma.asym0.prop90.Rda
limma_prad_asym0_80: $out_path/ext_analysis/prad.data.limma.asym0.prop80.Rda
limma_prad_asym0_70: $out_path/ext_analysis/prad.data.limma.asym0.prop70.Rda

limma_prad_asym1_90: $out_path/ext_analysis/prad.data.limma.asym1.prop90.Rda
limma_prad_asym1_80: $out_path/ext_analysis/prad.data.limma.asym1.prop80.Rda
limma_prad_asym1_70: $out_path/ext_analysis/prad.data.limma.asym1.prop70.Rda

limma_thca_asym0_90: $out_path/ext_analysis/thca.data.limma.asym0.prop90.Rda
limma_thca_asym0_80: $out_path/ext_analysis/thca.data.limma.asym0.prop80.Rda
limma_thca_asym0_70: $out_path/ext_analysis/thca.data.limma.asym0.prop70.Rda

limma_thca_asym1_90: $out_path/ext_analysis/thca.data.limma.asym1.prop90.Rda
limma_thca_asym1_80: $out_path/ext_analysis/thca.data.limma.asym1.prop80.Rda
limma_thca_asym1_70: $out_path/ext_analysis/thca.data.limma.asym1.prop70.Rda

#16S datasets

cdiff_edgeR: $out_path/ext_analysis/cdiff.data.edgeR.Rda

cdiff_limma: $out_path/ext_analysis/cdiff.data.limma.Rda
hiv_limma: $out_path/ext_analysis/hiv.data.limma.Rda
fre_limma: $out_path/ext_analysis/fre.data.limma.Rda
mar_limma: $out_path/ext_analysis/mar.data.limma.Rda






#rules to generate the output files
### DESeq
$out_path/ext_analysis/immuno.data.deseq.Rda : code/immuno_deseq.R code/des.fun.R
	Rscript 'code/immuno_deseq.R' 'code/des.fun.R'

$out_path/ext_analysis/brca.data.deseq.Rda : code/brca_deseq.R code/des.fun.R
	Rscript 'code/brca_deseq.R' 'code/des.fun.R'
	
$out_path/ext_analysis/kirc.data.deseq.Rda : code/kirc_deseq.R code/des.fun.R
	Rscript 'code/kirc_deseq.R' 'code/des.fun.R'

$out_path/ext_analysis/lihc.data.deseq.Rda : code/lihc_deseq.R code/des.fun.R
	Rscript 'code/lihc_deseq.R' 'code/des.fun.R'

$out_path/ext_analysis/luad.data.deseq.Rda : code/luad_deseq.R code/des.fun.R
	Rscript 'code/luad_deseq.R' 'code/des.fun.R'

$out_path/ext_analysis/prad.data.deseq.Rda : code/prad_deseq.R code/des.fun.R
	Rscript 'code/prad_deseq.R' 'code/des.fun.R'

$out_path/ext_analysis/thca.data.deseq.Rda : code/thca_deseq.R code/des.fun.R
	Rscript 'code/thca_deseq.R' 'code/des.fun.R'


### edgeR
$out_path/ext_analysis/immuno.data.edger.Rda : code/immuno_edgeR.R code/edg.fun.R
	Rscript 'code/immuno_edgeR.R' 'code/edg.fun.R'
	
$out_path/ext_analysis/brca.data.edger.Rda : code/brca_edgeR.R code/edg.fun.R
	Rscript 'code/brca_edgeR.R' 'code/edg.fun.R'
	
$out_path/ext_analysis/kirc.data.edger.Rda : code/kirc_edgeR.R code/edg.fun.R
	Rscript 'code/kirc_edgeR.R' 'code/edg.fun.R'

$out_path/ext_analysis/lihc.data.edger.Rda : code/lihc_edgeR.R code/edg.fun.R
	Rscript 'code/lihc_edgeR.R' 'code/edg.fun.R'
	
$out_path/ext_analysis/luad.data.edger.Rda : code/luad_edgeR.R code/edg.fun.R
	Rscript 'code/luad_edgeR.R' 'code/edg.fun.R'
	
$out_path/ext_analysis/prad.data.edger.Rda : code/prad_edgeR.R code/edg.fun.R
	Rscript 'code/prad_edgeR.R' 'code/edg.fun.R'

$out_path/ext_analysis/thca.data.edger.Rda : code/thca_edgeR.R code/edg.fun.R
	Rscript 'code/thca_edgeR.R' 'code/edg.fun.R'


### aldex2_0
$out_path/ext_analysis/immuno.data.aldex2_0.Rda : code/immuno_aldex2_0.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2_0.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/brca.data.aldex2_0.Rda : code/brca_aldex2_0.R code/ald2.fun.R
	Rscript 'code/brca_aldex2_0.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/kirc.data.aldex2_0.Rda : code/kirc_aldex2_0.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2_0.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/lihc.data.aldex2_0.Rda : code/lihc_aldex2_0.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2_0.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/luad.data.aldex2_0.Rda : code/luad_aldex2_0.R 
	Rscript 'code/luad_aldex2_0.R' 
	
$out_path/ext_analysis/prad.data.aldex2_0.Rda : code/prad_aldex2_0.R code/ald2.fun.R
	Rscript 'code/prad_aldex2_0.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/thca.data.aldex2_0.Rda : code/thca_aldex2_0.R code/ald2.fun.R
	Rscript 'code/thca_aldex2_0.R' 'code/ald2.fun.R'
	
### aldex2_1
$out_path/ext_analysis/immuno.data.aldex2_1.Rda : code/immuno_aldex2_1.R
	Rscript 'code/immuno_aldex2_1.R' 
	
$out_path/ext_analysis/brca.data.aldex2_1.Rda : code/brca_aldex2_1.R 
	Rscript 'code/brca_aldex2_1.R' 
	
$out_path/ext_analysis/kirc.data.aldex2_1.Rda : code/kirc_aldex2_1.R 
	Rscript 'code/kirc_aldex2_1.R'
	
$out_path/ext_analysis/lihc.data.aldex2_1.Rda : code/lihc_aldex2_1.R 
	Rscript 'code/lihc_aldex2_1.R' 
	
$out_path/ext_analysis/luad.data.aldex2_1.Rda : code/luad_aldex2_1.R 
	Rscript 'code/luad_aldex2_1.R' 
	
$out_path/ext_analysis/prad.data.aldex2_1.Rda : code/prad_aldex2_1.R 
	Rscript 'code/prad_aldex2_1.R' 
	
$out_path/ext_analysis/thca.data.aldex2_1.Rda : code/thca_aldex2_1.R 
	Rscript 'code/thca_aldex2_1.R' 
	
### aldex2_2
$out_path/ext_analysis/immuno.data.aldex2_2.Rda : code/immuno_aldex2_2.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2_2.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/brca.data.aldex2_2.Rda : code/brca_aldex2_2.R code/ald2.fun.R
	Rscript 'code/brca_aldex2_2.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/kirc.data.aldex2_2.Rda : code/kirc_aldex2_2.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2_2.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/lihc.data.aldex2_2.Rda : code/lihc_aldex2_2.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2_2.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/luad.data.aldex2_2.Rda : code/luad_aldex2_2.R code/ald2.fun.R
	Rscript 'code/luad_aldex2_2.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/prad.data.aldex2_2.Rda : code/prad_aldex2_2.R code/ald2.fun.R
	Rscript 'code/prad_aldex2_2.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/thca.data.aldex2_2.Rda : code/thca_aldex2_2.R code/ald2.fun.R
	Rscript 'code/thca_aldex2_2.R' 'code/ald2.fun.R'
	
### aldex2_5
$out_path/ext_analysis/immuno.data.aldex2_5.Rda : code/immuno_aldex2_5.R code/ald2.fun.R
	Rscript 'code/immuno_aldex2_5.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/brca.data.aldex2_5.Rda : code/brca_aldex2_5.R code/ald2.fun.R
	Rscript 'code/brca_aldex2_5.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/kirc.data.aldex2_5.Rda : code/kirc_aldex2_5.R code/ald2.fun.R
	Rscript 'code/kirc_aldex2_5.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/lihc.data.aldex2_5.Rda : code/lihc_aldex2_5.R code/ald2.fun.R
	Rscript 'code/lihc_aldex2_5.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/luad.data.aldex2_5.Rda : code/luad_aldex2_5.R code/ald2.fun.R
	Rscript 'code/luad_aldex2_5.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/prad.data.aldex2_5.Rda : code/prad_aldex2_5.R code/ald2.fun.R
	Rscript 'code/prad_aldex2_5.R' 'code/ald2.fun.R'
	
$out_path/ext_analysis/thca.data.aldex2_5.Rda : code/thca_aldex2_5.R code/ald2.fun.R
	Rscript 'code/thca_aldex2_5.R' 'code/ald2.fun.R'
	

###limma 	
$out_path/ext_analysis/immuno.data.limma.Rda : code/immuno_limma.R code/lim.fun.R
	Rscript 'code/immuno_limma.R' 'code/lim.fun.R'
	
$out_path/ext_analysis/brca.data.limma.Rda : code/brca_limma.R code/lim.fun.R
	Rscript 'code/brca_limma.R' 'code/lim.fun.R'
	
$out_path/ext_analysis/kirc.data.limma.Rda : code/kirc_limma.R code/lim.fun.R
	Rscript 'code/kirc_limma.R' 'code/lim.fun.R'
	
$out_path/ext_analysis/lihc.data.limma.Rda : code/lihc_limma.R code/lim.fun.R
	Rscript 'code/lihc_limma.R' 'code/lim.fun.R'

$out_path/ext_analysis/luad.data.limma.Rda : code/luad_limma.R code/lim.fun.R
	Rscript 'code/luad_limma.R' 'code/lim.fun.R'
	
$out_path/ext_analysis/prad.data.limma.Rda : code/prad_limma.R code/lim.fun.R
	Rscript 'code/prad_limma.R' 'code/lim.fun.R'
	
$out_path/ext_analysis/thca.data.limma.Rda : code/thca_limma.R code/lim.fun.R
	Rscript 'code/thca_limma.R' 'code/lim.fun.R'

#asymmetrical limma files
$out_path/ext_analysis/immuno.data.limma.asym1.Rda : code/immuno_limma_asym1.R
	Rscript 'code/immuno_limma_asym1.R'
	
$out_path/ext_analysis/immuno.data.limma.asym2.Rda : code/immuno_limma_asym2.R
	Rscript 'code/immuno_limma_asym2.R'
	
$out_path/ext_analysis/immuno.data.limma.asym3.Rda : code/immuno_limma_asym3.R
	Rscript 'code/immuno_limma_asym3.R'
	

$out_path/ext_analysis/brca.data.limma.asym1.Rda : code/brca_limma_asym1.R
	Rscript 'code/brca_limma_asym1.R'
	
$out_path/ext_analysis/brca.data.limma.asym2.Rda : code/brca_limma_asym2.R
	Rscript 'code/brca_limma_asym2.R'
	
$out_path/ext_analysis/brca.data.limma.asym3.Rda : code/brca_limma_asym3.R
	Rscript 'code/brca_limma_asym3.R'

$out_path/ext_analysis/kirc.data.limma.asym1.Rda : code/kirc_limma_asym1.R
	Rscript 'code/kirc_limma_asym1.R'
	
$out_path/ext_analysis/kirc.data.limma.asym2.Rda : code/kirc_limma_asym2.R
	Rscript 'code/kirc_limma_asym2.R'
	
$out_path/ext_analysis/kirc.data.limma.asym3.Rda : code/kirc_limma_asym3.R
	Rscript 'code/kirc_limma_asym3.R'

$out_path/ext_analysis/lihc.data.limma.asym1.Rda : code/lihc_limma_asym1.R
	Rscript 'code/lihc_limma_asym1.R'
	
$out_path/ext_analysis/lihc.data.limma.asym2.Rda : code/lihc_limma_asym2.R
	Rscript 'code/lihc_limma_asym2.R'
	
$out_path/ext_analysis/lihc.data.limma.asym3.Rda : code/lihc_limma_asym3.R
	Rscript 'code/lihc_limma_asym3.R'

$out_path/ext_analysis/luad.data.limma.asym1.Rda : code/luad_limma_asym1.R
	Rscript 'code/luad_limma_asym1.R'
	
$out_path/ext_analysis/luad.data.limma.asym2.Rda : code/luad_limma_asym2.R
	Rscript 'code/luad_limma_asym2.R'
	
$out_path/ext_analysis/luad.data.limma.asym3.Rda : code/luad_limma_asym3.R
	Rscript 'code/luad_limma_asym3.R'

$out_path/ext_analysis/prad.data.limma.asym1.Rda : code/prad_limma_asym1.R
	Rscript 'code/prad_limma_asym1.R'
	
$out_path/ext_analysis/prad.data.limma.asym2.Rda : code/prad_limma_asym2.R
	Rscript 'code/prad_limma_asym2.R'
	
$out_path/ext_analysis/prad.data.limma.asym3.Rda : code/prad_limma_asym3.R
	Rscript 'code/prad_limma_asym3.R'

$out_path/ext_analysis/thca.data.limma.asym1.Rda : code/thca_limma_asym1.R
	Rscript 'code/thca_limma_asym1.R'
	
$out_path/ext_analysis/thca.data.limma.asym2.Rda : code/thca_limma_asym2.R
	Rscript 'code/thca_limma_asym2.R'
	
$out_path/ext_analysis/thca.data.limma.asym3.Rda : code/thca_limma_asym3.R
	Rscript 'code/thca_limma_asym3.R'
	

### aldex3_0
$out_path/ext_analysis/immuno.data.aldex3_0.Rda : code/immuno_aldex3_0.R 
	Rscript 'code/immuno_aldex3_0.R' 
	
$out_path/ext_analysis/brca.data.aldex3_0.Rda : code/brca_aldex3_0.R 
	Rscript 'code/brca_aldex3_0.R' 
	
$out_path/ext_analysis/kirc.data.aldex3_0.Rda : code/kirc_aldex3_0.R 
	Rscript 'code/kirc_aldex3_0.R'
	
$out_path/ext_analysis/lihc.data.aldex3_0.Rda : code/lihc_aldex3_0.R 
	Rscript 'code/lihc_aldex3_0.R' 
	
$out_path/ext_analysis/luad.data.aldex3_0.Rda : code/luad_aldex3_0.R 
	Rscript 'code/luad_aldex3_0.R' 
	
$out_path/ext_analysis/prad.data.aldex3_0.Rda : code/prad_aldex3_0.R 
	Rscript 'code/prad_aldex3_0.R' 
	
$out_path/ext_analysis/thca.data.aldex3_0.Rda : code/thca_aldex3_0.R 
	Rscript 'code/thca_aldex3_0.R' 
	
### aldex3_1
$out_path/ext_analysis/immuno.data.aldex3_1.Rda : code/immuno_aldex3_1.R
	Rscript 'code/immuno_aldex3_1.R' 
	
$out_path/ext_analysis/brca.data.aldex3_1.Rda : code/brca_aldex3_1.R 
	Rscript 'code/brca_aldex3_1.R' 
	
$out_path/ext_analysis/kirc.data.aldex3_1.Rda : code/kirc_aldex3_1.R 
	Rscript 'code/kirc_aldex3_1.R'
	
$out_path/ext_analysis/lihc.data.aldex3_1.Rda : code/lihc_aldex3_1.R 
	Rscript 'code/lihc_aldex3_1.R' 
	
$out_path/ext_analysis/luad.data.aldex3_1.Rda : code/luad_aldex3_1.R 
	Rscript 'code/luad_aldex3_1.R' 
	
$out_path/ext_analysis/prad.data.aldex3_1.Rda : code/prad_aldex3_1.R 
	Rscript 'code/prad_aldex3_1.R' 
	
$out_path/ext_analysis/thca.data.aldex3_1.Rda : code/thca_aldex3_1.R 
	Rscript 'code/thca_aldex3_1.R' 
	
### aldex3_2
$out_path/ext_analysis/immuno.data.aldex3_2.Rda : code/immuno_aldex3_2.R 
	Rscript 'code/immuno_aldex3_2.R' 
	
$out_path/ext_analysis/brca.data.aldex3_2.Rda : code/brca_aldex3_2.R 
	Rscript 'code/brca_aldex3_2.R' 
	
$out_path/ext_analysis/kirc.data.aldex3_2.Rda : code/kirc_aldex3_2.R 
	Rscript 'code/kirc_aldex3_2.R' 
	
$out_path/ext_analysis/lihc.data.aldex3_2.Rda : code/lihc_aldex3_2.R 
	Rscript 'code/lihc_aldex3_2.R' 
	
$out_path/ext_analysis/luad.data.aldex3_2.Rda : code/luad_aldex3_2.R 
	Rscript 'code/luad_aldex3_2.R'
	
$out_path/ext_analysis/prad.data.aldex3_2.Rda : code/prad_aldex3_2.R 
	Rscript 'code/prad_aldex3_2.R' 
	
$out_path/ext_analysis/thca.data.aldex3_2.Rda : code/thca_aldex3_2.R 
	Rscript 'code/thca_aldex3_2.R'
	
### aldex3_5
$out_path/ext_analysis/immuno.data.aldex3_5.Rda : code/immuno_aldex3_5.R 
	Rscript 'code/immuno_aldex3_5.R' 
	
$out_path/ext_analysis/brca.data.aldex3_5.Rda : code/brca_aldex3_5.R 
	Rscript 'code/brca_aldex3_5.R' 

$out_path/ext_analysis/kirc.data.aldex3_5.Rda : code/kirc_aldex3_5.R 
	Rscript 'code/kirc_aldex3_5.R' 
	
$out_path/ext_analysis/lihc.data.aldex3_5.Rda : code/lihc_aldex3_5.R 
	Rscript 'code/lihc_aldex3_5.R' 
	
$out_path/ext_analysis/luad.data.aldex3_5.Rda : code/luad_aldex3_5.R 
	Rscript 'code/luad_aldex3_5.R' 
	
$out_path/ext_analysis/prad.data.aldex3_5.Rda : code/prad_aldex3_5.R 
	Rscript 'code/prad_aldex3_5.R' 
	
$out_path/ext_analysis/thca.data.aldex3_5.Rda : code/thca_aldex3_5.R 
	Rscript 'code/thca_aldex3_5.R'
	
#asymmetrical limma with different prop_null values

$out_path/ext_analysis/immuno.data.limma.asym0.prop90.Rda : code/immuno_limma_asym0_90.R
	Rscript 'code/immuno_limma_asym0_90.R'

$out_path/ext_analysis/immuno.data.limma.asym0.prop80.Rda : code/immuno_limma_asym0_80.R
	Rscript 'code/immuno_limma_asym0_80.R'

$out_path/ext_analysis/immuno.data.limma.asym0.prop70.Rda : code/immuno_limma_asym0_70.R
	Rscript 'code/immuno_limma_asym0_70.R'

$out_path/ext_analysis/immuno.data.limma.asym1.prop90.Rda : code/immuno_limma_asym1_90.R
	Rscript 'code/immuno_limma_asym1_90.R'

$out_path/ext_analysis/immuno.data.limma.asym1.prop80.Rda : code/immuno_limma_asym1_80.R
	Rscript 'code/immuno_limma_asym1_80.R'

$out_path/ext_analysis/immuno.data.limma.asym1.prop70.Rda : code/immuno_limma_asym1_70.R
	Rscript 'code/immuno_limma_asym1_70.R'

$out_path/ext_analysis/brca.data.limma.asym0.prop90.Rda : code/brca_limma_asym0_90.R
	Rscript 'code/brca_limma_asym0_90.R'

$out_path/ext_analysis/brca.data.limma.asym0.prop80.Rda : code/brca_limma_asym0_80.R
	Rscript 'code/brca_limma_asym0_80.R'

$out_path/ext_analysis/brca.data.limma.asym0.prop70.Rda : code/brca_limma_asym0_70.R
	Rscript 'code/brca_limma_asym0_70.R'

$out_path/ext_analysis/brca.data.limma.asym1.prop90.Rda : code/brca_limma_asym1_90.R
	Rscript 'code/brca_limma_asym1_90.R'

$out_path/ext_analysis/brca.data.limma.asym1.prop80.Rda : code/brca_limma_asym1_80.R
	Rscript 'code/brca_limma_asym1_80.R'

$out_path/ext_analysis/brca.data.limma.asym1.prop70.Rda : code/brca_limma_asym1_70.R
	Rscript 'code/brca_limma_asym1_70.R'
	
$out_path/ext_analysis/kirc.data.limma.asym0.prop90.Rda : code/kirc_limma_asym0_90.R
	Rscript 'code/kirc_limma_asym0_90.R'

$out_path/ext_analysis/kirc.data.limma.asym0.prop80.Rda : code/kirc_limma_asym0_80.R
	Rscript 'code/kirc_limma_asym0_80.R'

$out_path/ext_analysis/kirc.data.limma.asym0.prop70.Rda : code/kirc_limma_asym0_70.R
	Rscript 'code/kirc_limma_asym0_70.R'

$out_path/ext_analysis/kirc.data.limma.asym1.prop90.Rda : code/kirc_limma_asym1_90.R
	Rscript 'code/kirc_limma_asym1_90.R'

$out_path/ext_analysis/kirc.data.limma.asym1.prop80.Rda : code/kirc_limma_asym1_80.R
	Rscript 'code/kirc_limma_asym1_80.R'

$out_path/ext_analysis/kirc.data.limma.asym1.prop70.Rda : code/kirc_limma_asym1_70.R
	Rscript 'code/kirc_limma_asym1_70.R'

$out_path/ext_analysis/lihc.data.limma.asym0.prop90.Rda : code/lihc_limma_asym0_90.R
	Rscript 'code/lihc_limma_asym0_90.R'

$out_path/ext_analysis/lihc.data.limma.asym0.prop80.Rda : code/lihc_limma_asym0_80.R
	Rscript 'code/lihc_limma_asym0_80.R'

$out_path/ext_analysis/lihc.data.limma.asym0.prop70.Rda : code/lihc_limma_asym0_70.R
	Rscript 'code/lihc_limma_asym0_70.R'

$out_path/ext_analysis/lihc.data.limma.asym1.prop90.Rda : code/lihc_limma_asym1_90.R
	Rscript 'code/lihc_limma_asym1_90.R'

$out_path/ext_analysis/lihc.data.limma.asym1.prop80.Rda : code/lihc_limma_asym1_80.R
	Rscript 'code/lihc_limma_asym1_80.R'

$out_path/ext_analysis/lihc.data.limma.asym1.prop70.Rda : code/lihc_limma_asym1_70.R
	Rscript 'code/lihc_limma_asym1_70.R'
	
$out_path/ext_analysis/luad.data.limma.asym0.prop90.Rda : code/luad_limma_asym0_90.R
	Rscript 'code/luad_limma_asym0_90.R'

$out_path/ext_analysis/luad.data.limma.asym0.prop80.Rda : code/luad_limma_asym0_80.R
	Rscript 'code/luad_limma_asym0_80.R'

$out_path/ext_analysis/luad.data.limma.asym0.prop70.Rda : code/luad_limma_asym0_70.R
	Rscript 'code/luad_limma_asym0_70.R'

$out_path/ext_analysis/luad.data.limma.asym1.prop90.Rda : code/luad_limma_asym1_90.R
	Rscript 'code/luad_limma_asym1_90.R'

$out_path/ext_analysis/luad.data.limma.asym1.prop80.Rda : code/luad_limma_asym1_80.R
	Rscript 'code/luad_limma_asym1_80.R'

$out_path/ext_analysis/luad.data.limma.asym1.prop70.Rda : code/luad_limma_asym1_70.R
	Rscript 'code/luad_limma_asym1_70.R'

$out_path/ext_analysis/prad.data.limma.asym0.prop90.Rda : code/prad_limma_asym0_90.R
	Rscript 'code/prad_limma_asym0_90.R'

$out_path/ext_analysis/prad.data.limma.asym0.prop80.Rda : code/prad_limma_asym0_80.R
	Rscript 'code/prad_limma_asym0_80.R'

$out_path/ext_analysis/prad.data.limma.asym0.prop70.Rda : code/prad_limma_asym0_70.R
	Rscript 'code/prad_limma_asym0_70.R'

$out_path/ext_analysis/prad.data.limma.asym1.prop90.Rda : code/prad_limma_asym1_90.R
	Rscript 'code/prad_limma_asym1_90.R'

$out_path/ext_analysis/prad.data.limma.asym1.prop80.Rda : code/prad_limma_asym1_80.R
	Rscript 'code/prad_limma_asym1_80.R'

$out_path/ext_analysis/prad.data.limma.asym1.prop70.Rda : code/prad_limma_asym1_70.R
	Rscript 'code/prad_limma_asym1_70.R'
	
$out_path/ext_analysis/thca.data.limma.asym0.prop90.Rda : code/thca_limma_asym0_90.R
	Rscript 'code/thca_limma_asym0_90.R'

$out_path/ext_analysis/thca.data.limma.asym0.prop80.Rda : code/thca_limma_asym0_80.R
	Rscript 'code/thca_limma_asym0_80.R'

$out_path/ext_analysis/thca.data.limma.asym0.prop70.Rda : code/thca_limma_asym0_70.R
	Rscript 'code/prad_limma_asym0_70.R'

$out_path/ext_analysis/thca.data.limma.asym1.prop90.Rda : code/thca_limma_asym1_90.R
	Rscript 'code/thca_limma_asym1_90.R'

$out_path/ext_analysis/thca.data.limma.asym1.prop80.Rda : code/thca_limma_asym1_80.R
	Rscript 'code/thca_limma_asym1_80.R'

$out_path/ext_analysis/thca.data.limma.asym1.prop70.Rda : code/thca_limma_asym1_70.R
	Rscript 'code/thca_limma_asym1_70.R'
	
#16S datasets

$out_path/ext_analysis/cdiff.data.edgeR.Rda : code/cdiff_edgeR.R
	Rscript 'code/cdiff_edgeR.R'

$out_path/ext_analysis/cdiff.data.limma.Rda : code/cdiff_limma.R
	Rscript 'code/cdiff_limma.R'
	
$out_path/ext_analysis/fre.data.limma.Rda : code/freshwater_limma.R
	Rscript 'code/freshwater_limma.R'
	
$out_path/ext_analysis/hiv.data.limma.Rda : code/humanhiv_limma.R
	Rscript 'code/humanhiv_limma.R'
	
$out_path/ext_analysis/mar.data.limma.Rda : code/marinesed_limma.R
	Rscript 'code/marinesed_limma.R'

#DESeq clean files
clean_immuno_DESeq:
	rm $out_path/ext_analysis/immuno.data.deseq.Rda
clean_brca_DESeq:
	rm $out_path/ext_analysis/brca.data.deseq.Rda
clean_kirc_DESeq:
	rm $out_path/ext_analysis/kirc.data.deseq.Rda
clean_lihc_DESeq:
	rm $out_path/ext_analysis/lihc.data.deseq.Rda
clean_luad_DESeq:
	rm $out_path/ext_analysis/luad.data.deseq.Rda
clean_prad_DESeq:
	rm $out_path/ext_analysis/prad.data.deseq.Rda
clean_thca_DESeq:
	rm $out_path/ext_analysis/thca.data.deseq.Rda

#edgeR clean files
clean_immuno_edgeR:
	rm $out_path/ext_analysis/immuno.data.edger.Rda
clean_brca_edgeR:
	rm $out_path/ext_analysis/brca.data.edger.Rda
clean_kirc_edgeR:
	rm $out_path/ext_analysis/kirc.data.edger.Rda
clean_lihc_edgeR:
	rm $out_path/ext_analysis/lihc.data.edger.Rda
clean_luad_edgeR:
	rm $out_path/ext_analysis/luad.data.edger.Rda
clean_prad_edgeR:
	rm $out_path/ext_analysis/prad.data.edger.Rda
clean_thca_edgeR:
	rm $out_path/ext_analysis/thca.data.edger.Rda
	
#aldex clean files
clean_immuno_aldex2_0:
	rm $out_path/ext_analysis/immuno.data.aldex2_0.Rda
clean_immuno_aldex2_1:
	rm $out_path/ext_analysis/immuno.data.aldex2_1.Rda
clean_immuno_aldex2_2:
	rm $out_path/ext_analysis/immuno.data.aldex2_2.Rda
clean_immuno_aldex2_5:
	rm $out_path/ext_analysis/immuno.data.aldex2_5.Rda
clean_brca_aldex2_0:
	rm $out_path/ext_analysis/brca.data.aldex2_0.Rda
clean_brca_aldex2_1:
	rm $out_path/ext_analysis/brca.data.aldex2_1.Rda
clean_brca_aldex2_2:
	rm $out_path/ext_analysis/brca.data.aldex2_2.Rda
clean_brca_aldex2_5:
	rm $out_path/ext_analysis/brca.data.aldex2_5.Rda
clean_kirc_aldex2_0:
	rm $out_path/ext_analysis/kirc.data.aldex2_0.Rda
clean_kirc_aldex2_1:
	rm $out_path/ext_analysis/kirc.data.aldex2_1.Rda
clean_kirc_aldex2_2:
	rm $out_path/ext_analysis/kirc.data.aldex2_2.Rda
clean_kirc_aldex2_5:
	rm $out_path/ext_analysis/kirc.data.aldex2_5.Rda
clean_lihc_aldex2_0:
	rm $out_path/ext_analysis/lihc.data.aldex2_0.Rda
clean_lihc_aldex2_1:
	rm $out_path/ext_analysis/lihc.data.aldex2_1.Rda
clean_lihc_aldex2_2:
	rm $out_path/ext_analysis/lihc.data.aldex2_2.Rda
clean_lihc_aldex2_5:
	rm $out_path/ext_analysis/lihc.data.aldex2_5.Rda
clean_luad_aldex2_0:
	rm $out_path/ext_analysis/luad.data.aldex2_0.Rda
clean_luad_aldex2_1:
	rm $out_path/ext_analysis/luad.data.aldex2_1.Rda
clean_luad_aldex2_2:
	rm $out_path/ext_analysis/luad.data.aldex2_2.Rda
clean_luad_aldex2_5:
	rm $out_path/ext_analysis/luad.data.aldex2_5.Rda
clean_prad_aldex2_0:
	rm $out_path/ext_analysis/prad.data.aldex2_0.Rda
clean_prad_aldex2_1:
	rm $out_path/ext_analysis/prad.data.aldex2_1.Rda
clean_prad_aldex2_2:
	rm $out_path/ext_analysis/prad.data.aldex2_2.Rda
clean_prad_aldex2_5:
	rm $out_path/ext_analysis/prad.data.aldex2_5.Rda
clean_thca_aldex2_0:
	rm $out_path/ext_analysis/thca.data.aldex2_0.Rda
clean_thca_aldex2_1:
	rm $out_path/ext_analysis/thca.data.aldex2_1.Rda
clean_thca_aldex2_2:
	rm $out_path/ext_analysis/thca.data.aldex2_2.Rda
clean_thca_aldex2_5:
	rm $out_path/ext_analysis/thca.data.aldex2_5.Rda

#limma clean files
clean_immuno_limma:
	rm $out_path/ext_analysis/immuno.data.limma.Rda
clean_brca_limma:
	rm $out_path/ext_analysis/brca.data.limma.Rda
clean_kirc_limma:
	rm $out_path/ext_analysis/kirc.data.limma.Rda
clean_lihc_limma:
	rm $out_path/ext_analysis/lihc.data.limma.Rda
clean_luad_limma:
	rm $out_path/ext_analysis/luad.data.limma.Rda
clean_prad_limma:
	rm $out_path/ext_analysis/prad.data.limma.Rda
clean_thca_limma:
	rm $out_path/ext_analysis/thca.data.limma.Rda

	
	
