#makes the file .Rda if make all is called
#all: analysis/test2.Rda

# define output directory
out_path="/Volumes/data2/andreea"



# define final target objects for each tool & dataset - no asymmetry

### DESeq2
DESeq_brca: $out_path/ext_analysis/brca.data.deseq.Rda
DESeq_cdiff: $out_path/ext_analysis/cdiff.data.deseq.Rda
DESeq_fre: $out_path/ext_analysis/fre.data.deseq.Rda
DESeq_hiv: $out_path/ext_analysis/hiv.data.deseq.Rda
DESeq_immuno: $out_path/ext_analysis/immuno.data.deseq.Rda
DESeq_kirc: $out_path/ext_analysis/kirc.data.deseq.Rda
DESeq_lihc: $out_path/ext_analysis/lihc.data.deseq.Rda
DESeq_luad: $out_path/ext_analysis/luad.data.deseq.Rda
DESeq_mar: $out_path/ext_analysis/mar.data.deseq.Rda
DESeq_mts: $out_path/ext_analysis/mts.data.deseq.Rda
DESeq_oral: $out_path/ext_analysis/oral.data.deseq.Rda
DESeq_prad: $out_path/ext_analysis/prad.data.deseq.Rda
DESeq_sc: $out_path/ext_analysis/sccyto.data.deseq.Rda
DESeq_thca: $out_path/ext_analysis/thca.data.deseq.Rda

### edgeR
edgeR_brca: $out_path/ext_analysis/brca.data.edger.Rda
edgeR_cdiff: $out_path/ext_analysis/cdiff.data.edger.Rda
edgeR_fre: $out_path/ext_analysis/fre.data.edger.Rda
edgeR_hiv: $out_path/ext_analysis/hiv.data.edger.Rda
edgeR_immuno: $out_path/ext_analysis/immuno.data.edger.Rda
edgeR_kirc: $out_path/ext_analysis/kirc.data.edger.Rda
edgeR_lihc: $out_path/ext_analysis/lihc.data.edger.Rda
edgeR_luad: $out_path/ext_analysis/luad.data.edger.Rda
edgeR_mar: $out_path/ext_analysis/mar.data.edger.Rda
edgeR_mts: $out_path/ext_analysis/mts.data.edger.Rda
edgeR_oral: $out_path/ext_analysis/oral.data.edger.Rda
edgeR_prad: $out_path/ext_analysis/prad.data.edger.Rda
edgeR_sc: $out_path/ext_analysis/sccyto.data.edger.Rda
edgeR_thca: $out_path/ext_analysis/thca.data.edger.Rda

### Limma
limma_brca: $out_path/ext_analysis/brca.data.limma.Rda
limma_cdiff: $out_path/ext_analysis/cdiff.data.limma.Rda
limma_fre: $out_path/ext_analysis/fre.data.limma.Rda
limma_hiv: $out_path/ext_analysis/hiv.data.limma.Rda
limma_immuno: $out_path/ext_analysis/immuno.data.limma.Rda
limma_kirc: $out_path/ext_analysis/kirc.data.limma.Rda
limma_lihc: $out_path/ext_analysis/lihc.data.limma.Rda
limma_luad: $out_path/ext_analysis/luad.data.limma.Rda
limma_mar: $out_path/ext_analysis/mar.data.limma.Rda
limma_mts: $out_path/ext_analysis/mts.data.limma.Rda
limma_oral: $out_path/ext_analysis/oral.data.limma.Rda
limma_prad: $out_path/ext_analysis/prad.data.limma.Rda
limma_sccyto: $out_path/ext_analysis/sccyto.data.limma.Rda
limma_thca: $out_path/ext_analysis/thca.data.limma.Rda

### ALDEx2 gamma = 1e-3
ALDEx2_0_brca: $out_path/ext_analysis/brca.data.aldex2_0.Rda
ALDEx2_0_cdiff: $out_path/ext_analysis/cdiff.data.aldex2_0.Rda
ALDEx2_0_fre: $out_path/ext_analysis/fre.data.aldex2_0.Rda
ALDEx2_0_hiv: $out_path/ext_analysis/hiv.data.aldex2_0.Rda
ALDEx2_0_immuno: $out_path/ext_analysis/immuno.data.aldex2_0.Rda
ALDEx2_0_kirc: $out_path/ext_analysis/kirc.data.aldex2_0.Rda
ALDEx2_0_lihc: $out_path/ext_analysis/lihc.data.aldex2_0.Rda
ALDEx2_0_luad: $out_path/ext_analysis/luad.data.aldex2_0.Rda
ALDEx2_0_mar: $out_path/ext_analysis/mar.data.aldex2_0.Rda
ALDEx2_0_mts: $out_path/ext_analysis/mts.data.aldex2_0.Rda
ALDEx2_0_oral: $out_path/ext_analysis/oral.data.aldex2_0.Rda
ALDEx2_0_prad: $out_path/ext_analysis/prad.data.aldex2_0.Rda
ALDEx2_0_sccyto: $out_path/ext_analysis/sccyto.data.aldex2_0.Rda
ALDEx2_0_thca: $out_path/ext_analysis/thca.data.aldex2_0.Rda

### ALDEx2 gamma = 0.1
ALDEx2_1_brca: $out_path/ext_analysis/brca.data.aldex2_1.Rda
ALDEx2_1_cdiff: $out_path/ext_analysis/cdiff.data.aldex2_1.Rda
ALDEx2_1_fre: $out_path/ext_analysis/fre.data.aldex2_1.Rda
ALDEx2_1_hiv: $out_path/ext_analysis/hiv.data.aldex2_1.Rda
ALDEx2_1_immuno: $out_path/ext_analysis/immuno.data.aldex2_1.Rda
ALDEx2_1_kirc: $out_path/ext_analysis/kirc.data.aldex2_1.Rda
ALDEx2_1_lihc: $out_path/ext_analysis/lihc.data.aldex2_1.Rda
ALDEx2_1_luad: $out_path/ext_analysis/luad.data.aldex2_1.Rda
ALDEx2_1_mar: $out_path/ext_analysis/mar.data.aldex2_1.Rda
ALDEx2_1_mts: $out_path/ext_analysis/mts.data.aldex2_1.Rda
ALDEx2_1_oral: $out_path/ext_analysis/oral.data.aldex2_1.Rda
ALDEx2_1_prad: $out_path/ext_analysis/prad.data.aldex2_1.Rda
ALDEx2_1_sccyto: $out_path/ext_analysis/sccyto.data.aldex2_1.Rda
ALDEx2_1_thca: $out_path/ext_analysis/thca.data.aldex2_1.Rda

### ALDEx2 gamma = 0.2
ALDEx2_2_brca: $out_path/ext_analysis/brca.data.aldex2_2.Rda
ALDEx2_2_cdiff: $out_path/ext_analysis/cdiff.data.aldex2_2.Rda
ALDEx2_2_fre: $out_path/ext_analysis/fre.data.aldex2_2.Rda
ALDEx2_2_hiv: $out_path/ext_analysis/hiv.data.aldex2_2.Rda
ALDEx2_2_immuno: $out_path/ext_analysis/immuno.data.aldex2_2.Rda
ALDEx2_2_kirc: $out_path/ext_analysis/kirc.data.aldex2_2.Rda
ALDEx2_2_lihc: $out_path/ext_analysis/lihc.data.aldex2_2.Rda
ALDEx2_2_luad: $out_path/ext_analysis/luad.data.aldex2_2.Rda
ALDEx2_2_mar: $out_path/ext_analysis/mar.data.aldex2_2.Rda
ALDEx2_2_mts: $out_path/ext_analysis/mts.data.aldex2_2.Rda
ALDEx2_2_oral: $out_path/ext_analysis/oral.data.aldex2_2.Rda
ALDEx2_2_prad: $out_path/ext_analysis/prad.data.aldex2_2.Rda
ALDEx2_2_sccyto: $out_path/ext_analysis/sccyto.data.aldex2_2.Rda
ALDEx2_2_thca: $out_path/ext_analysis/thca.data.aldex2_2.Rda

### ALDEx2 gamma = 0.3
ALDEx2_3_brca: $out_path/ext_analysis/brca.data.aldex2_3.Rda
ALDEx2_3_cdiff: $out_path/ext_analysis/cdiff.data.aldex2_3.Rda
ALDEx2_3_fre: $out_path/ext_analysis/fre.data.aldex2_3.Rda
ALDEx2_3_hiv: $out_path/ext_analysis/hiv.data.aldex2_3.Rda
ALDEx2_3_immuno: $out_path/ext_analysis/immuno.data.aldex2_3.Rda
ALDEx2_3_kirc: $out_path/ext_analysis/kirc.data.aldex2_3.Rda
ALDEx2_3_lihc: $out_path/ext_analysis/lihc.data.aldex2_3.Rda
ALDEx2_3_luad: $out_path/ext_analysis/luad.data.aldex2_3.Rda
ALDEx2_3_mar: $out_path/ext_analysis/mar.data.aldex2_3.Rda
ALDEx2_3_mts: $out_path/ext_analysis/mts.data.aldex2_2.Rda
ALDEx2_3_oral: $out_path/ext_analysis/oral.data.aldex2_3.Rda
ALDEx2_3_prad: $out_path/ext_analysis/prad.data.aldex2_3.Rda
ALDEx2_3_sccyto: $out_path/ext_analysis/sccyto.data.aldex2_3.Rda
ALDEx2_3_thca: $out_path/ext_analysis/thca.data.aldex2_3.Rda

### ALDEx2 gamma = 0.4
ALDEx2_4_brca: $out_path/ext_analysis/brca.data.aldex2_4.Rda
ALDEx2_4_cdiff: $out_path/ext_analysis/cdiff.data.aldex2_4.Rda
ALDEx2_4_fre: $out_path/ext_analysis/fre.data.aldex2_4.Rda
ALDEx2_4_hiv: $out_path/ext_analysis/hiv.data.aldex2_4.Rda
ALDEx2_4_immuno: $out_path/ext_analysis/immuno.data.aldex2_4.Rda
ALDEx2_4_kirc: $out_path/ext_analysis/kirc.data.aldex2_4.Rda
ALDEx2_4_lihc: $out_path/ext_analysis/lihc.data.aldex2_4.Rda
ALDEx2_4_luad: $out_path/ext_analysis/luad.data.aldex2_4.Rda
ALDEx2_4_mar: $out_path/ext_analysis/mar.data.aldex2_4.Rda
ALDEx2_4_mts: $out_path/ext_analysis/mts.data.aldex2_2.Rda
ALDEx2_4_oral: $out_path/ext_analysis/oral.data.aldex2_4.Rda
ALDEx2_4_prad: $out_path/ext_analysis/prad.data.aldex2_4.Rda
ALDEx2_4_sccyto: $out_path/ext_analysis/sccyto.data.aldex2_4.Rda
ALDEx2_4_thca: $out_path/ext_analysis/thca.data.aldex2_4.Rda

### ALDEx2 gamma = 0.5
ALDEx2_5_brca: $out_path/ext_analysis/brca.data.aldex2_5.Rda
ALDEx2_5_cdiff: $out_path/ext_analysis/cdiff.data.aldex2_5.Rda
ALDEx2_5_fre: $out_path/ext_analysis/fre.data.aldex2_5.Rda
ALDEx2_5_hiv: $out_path/ext_analysis/hiv.data.aldex2_5.Rda
ALDEx2_5_immuno: $out_path/ext_analysis/immuno.data.aldex2_5.Rda
ALDEx2_5_kirc: $out_path/ext_analysis/kirc.data.aldex2_5.Rda
ALDEx2_5_lihc: $out_path/ext_analysis/lihc.data.aldex2_5.Rda
ALDEx2_5_luad: $out_path/ext_analysis/luad.data.aldex2_5.Rda
ALDEx2_5_mar: $out_path/ext_analysis/mar.data.aldex2_5.Rda
ALDEx2_5_mts: $out_path/ext_analysis/mts.data.aldex2_5.Rda
ALDEx2_5_oral: $out_path/ext_analysis/oral.data.aldex2_5.Rda
ALDEx2_5_prad: $out_path/ext_analysis/prad.data.aldex2_5.Rda
ALDEx2_5_sccyto: $out_path/ext_analysis/sccyto.data.aldex2_5.Rda
ALDEx2_5_thca: $out_path/ext_analysis/thca.data.aldex2_5.Rda

### ALDEx3 gamma = 1e-3
ALDEx3_0_brca: $out_path/ext_analysis/brca.data.aldex3_0.Rda
ALDEx3_0_cdiff: $out_path/ext_analysis/cdiff.data.aldex3_0.Rda
ALDEx3_0_fre: $out_path/ext_analysis/fre.data.aldex3_0.Rda
ALDEx3_0_hiv: $out_path/ext_analysis/hiv.data.aldex3_0.Rda
ALDEx3_0_immuno: $out_path/ext_analysis/immuno.data.aldex3_0.Rda
ALDEx3_0_kirc: $out_path/ext_analysis/kirc.data.aldex3_0.Rda
ALDEx3_0_lihc: $out_path/ext_analysis/lihc.data.aldex3_0.Rda
ALDEx3_0_luad: $out_path/ext_analysis/luad.data.aldex3_0.Rda
ALDEx3_0_mar: $out_path/ext_analysis/mar.data.aldex3_0.Rda
ALDEx3_0_mts: $out_path/ext_analysis/mts.data.aldex3_0.Rda
ALDEx3_0_oral: $out_path/ext_analysis/oral.data.aldex3_0.Rda
ALDEx3_0_prad: $out_path/ext_analysis/prad.data.aldex3_0.Rda
ALDEx3_0_sccyto: $out_path/ext_analysis/sccyto.data.aldex3_0.Rda
ALDEx3_0_thca: $out_path/ext_analysis/thca.data.aldex3_0.Rda

### ALDEx3 gamma = 0.1
ALDEx3_1_brca: $out_path/ext_analysis/brca.data.aldex3_1.Rda
ALDEx3_1_cdiff: $out_path/ext_analysis/cdiff.data.aldex3_1.Rda
ALDEx3_1_fre: $out_path/ext_analysis/fre.data.aldex3_1.Rda
ALDEx3_1_hiv: $out_path/ext_analysis/hiv.data.aldex3_1.Rda
ALDEx3_1_immuno: $out_path/ext_analysis/immuno.data.aldex3_1.Rda
ALDEx3_1_kirc: $out_path/ext_analysis/kirc.data.aldex3_1.Rda
ALDEx3_1_lihc: $out_path/ext_analysis/lihc.data.aldex3_1.Rda
ALDEx3_1_luad: $out_path/ext_analysis/luad.data.aldex3_1.Rda
ALDEx3_1_mar: $out_path/ext_analysis/mar.data.aldex3_1.Rda
ALDEx3_1_mts: $out_path/ext_analysis/mts.data.aldex3_1.Rda
ALDEx3_1_oral: $out_path/ext_analysis/oral.data.aldex3_1.Rda
ALDEx3_1_prad: $out_path/ext_analysis/prad.data.aldex3_1.Rda
ALDEx3_1_sccyto: $out_path/ext_analysis/sccyto.data.aldex3_1.Rda
ALDEx3_1_thca: $out_path/ext_analysis/thca.data.aldex3_1.Rda

### ALDEx3 gamma = 0.2
ALDEx3_2_brca: $out_path/ext_analysis/brca.data.aldex3_2.Rda
ALDEx3_2_cdiff: $out_path/ext_analysis/cdiff.data.aldex3_2.Rda
ALDEx3_2_fre: $out_path/ext_analysis/fre.data.aldex3_2.Rda
ALDEx3_2_hiv: $out_path/ext_analysis/hiv.data.aldex3_2.Rda
ALDEx3_2_immuno: $out_path/ext_analysis/immuno.data.aldex3_2.Rda
ALDEx3_2_kirc: $out_path/ext_analysis/kirc.data.aldex3_2.Rda
ALDEx3_2_lihc: $out_path/ext_analysis/lihc.data.aldex3_2.Rda
ALDEx3_2_luad: $out_path/ext_analysis/luad.data.aldex3_2.Rda
ALDEx3_2_mar: $out_path/ext_analysis/mar.data.aldex3_2.Rda
ALDEx3_2_mts: $out_path/ext_analysis/mts.data.aldex3_2.Rda
ALDEx3_2_oral: $out_path/ext_analysis/oral.data.aldex3_2.Rda
ALDEx3_2_prad: $out_path/ext_analysis/prad.data.aldex3_2.Rda
ALDEx3_2_sccyto: $out_path/ext_analysis/sccyto.data.aldex3_2.Rda
ALDEx3_2_thca: $out_path/ext_analysis/thca.data.aldex3_2.Rda

### ALDEx3 gamma = 0.3
ALDEx3_3_brca: $out_path/ext_analysis/brca.data.aldex3_3.Rda
ALDEx3_3_cdiff: $out_path/ext_analysis/cdiff.data.aldex3_3.Rda
ALDEx3_3_fre: $out_path/ext_analysis/fre.data.aldex3_3.Rda
ALDEx3_3_hiv: $out_path/ext_analysis/hiv.data.aldex3_3.Rda
ALDEx3_3_immuno: $out_path/ext_analysis/immuno.data.aldex3_3.Rda
ALDEx3_3_kirc: $out_path/ext_analysis/kirc.data.aldex3_3.Rda
ALDEx3_3_lihc: $out_path/ext_analysis/lihc.data.aldex3_3.Rda
ALDEx3_3_luad: $out_path/ext_analysis/luad.data.aldex3_3.Rda
ALDEx3_3_mar: $out_path/ext_analysis/mar.data.aldex3_3.Rda
ALDEx3_3_mts: $out_path/ext_analysis/mts.data.aldex3_3.Rda
ALDEx3_3_oral: $out_path/ext_analysis/oral.data.aldex3_3.Rda
ALDEx3_3_prad: $out_path/ext_analysis/prad.data.aldex3_3.Rda
ALDEx3_3_sccyto: $out_path/ext_analysis/sccyto.data.aldex3_3.Rda
ALDEx3_3_thca: $out_path/ext_analysis/thca.data.aldex3_3.Rda

### ALDEx3 gamma = 0.4
ALDEx3_4_brca: $out_path/ext_analysis/brca.data.aldex3_4.Rda
ALDEx3_4_cdiff: $out_path/ext_analysis/cdiff.data.aldex3_4.Rda
ALDEx3_4_fre: $out_path/ext_analysis/fre.data.aldex3_4.Rda
ALDEx3_4_hiv: $out_path/ext_analysis/hiv.data.aldex3_4.Rda
ALDEx3_4_immuno: $out_path/ext_analysis/immuno.data.aldex3_4.Rda
ALDEx3_4_kirc: $out_path/ext_analysis/kirc.data.aldex3_4.Rda
ALDEx3_4_lihc: $out_path/ext_analysis/lihc.data.aldex3_4.Rda
ALDEx3_4_luad: $out_path/ext_analysis/luad.data.aldex3_4.Rda
ALDEx3_4_mar: $out_path/ext_analysis/mar.data.aldex3_4.Rda
ALDEx3_4_mts: $out_path/ext_analysis/mts.data.aldex3_4.Rda
ALDEx3_4_oral: $out_path/ext_analysis/oral.data.aldex3_4.Rda
ALDEx3_4_prad: $out_path/ext_analysis/prad.data.aldex3_4.Rda
ALDEx3_4_sccyto: $out_path/ext_analysis/sccyto.data.aldex3_4.Rda
ALDEx3_4_thca: $out_path/ext_analysis/thca.data.aldex3_4.Rda

### ALDEx3 gamma = 0.5
ALDEx3_5_brca: $out_path/ext_analysis/brca.data.aldex3_5.Rda
ALDEx3_5_cdiff: $out_path/ext_analysis/cdiff.data.aldex3_5.Rda
ALDEx3_5_fre: $out_path/ext_analysis/fre.data.aldex3_5.Rda
ALDEx3_5_hiv: $out_path/ext_analysis/hiv.data.aldex3_5.Rda
ALDEx3_5_immuno: $out_path/ext_analysis/immuno.data.aldex3_5.Rda
ALDEx3_5_kirc: $out_path/ext_analysis/kirc.data.aldex3_5.Rda
ALDEx3_5_lihc: $out_path/ext_analysis/lihc.data.aldex3_5.Rda
ALDEx3_5_luad: $out_path/ext_analysis/luad.data.aldex3_5.Rda
ALDEx3_5_mar: $out_path/ext_analysis/mar.data.aldex3_5.Rda
ALDEx3_5_mts: $out_path/ext_analysis/mts.data.aldex3_5.Rda
ALDEx3_5_oral: $out_path/ext_analysis/oral.data.aldex3_5.Rda
ALDEx3_5_prad: $out_path/ext_analysis/prad.data.aldex3_5.Rda
ALDEx3_5_sccyto: $out_path/ext_analysis/sccyto.data.aldex3_5.Rda
ALDEx3_5_thca: $out_path/ext_analysis/thca.data.aldex3_5.Rda

### limma: asymmetry around mean of normal distribution (thinning)
limma_brca_asym1: $out_path/ext_analysis/brca.data.limma.asym1.Rda
limma_brca_asym2: $out_path/ext_analysis/brca.data.limma.asym2.Rda
limma_brca_asym3: $out_path/ext_analysis/brca.data.limma.asym3.Rda
limma_immuno_asym1: $out_path/ext_analysis/immuno.data.limma.asym1.Rda
limma_immuno_asym2: $out_path/ext_analysis/immuno.data.limma.asym2.Rda
limma_immuno_asym3: $out_path/ext_analysis/immuno.data.limma.asym3.Rda
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


### limma: asymmetry around mean of normal distribution and proportion of genes NOT altered (thinning)
limma_brca_asym0_90: $out_path/ext_analysis/brca.data.limma.asym0.prop90.Rda
limma_brca_asym0_80: $out_path/ext_analysis/brca.data.limma.asym0.prop80.Rda
limma_brca_asym0_70: $out_path/ext_analysis/brca.data.limma.asym0.prop70.Rda
limma_brca_asym1_90: $out_path/ext_analysis/brca.data.limma.asym1.prop90.Rda
limma_brca_asym1_80: $out_path/ext_analysis/brca.data.limma.asym1.prop80.Rda
limma_brca_asym1_70: $out_path/ext_analysis/brca.data.limma.asym1.prop70.Rda

limma_immuno_asym0_90: $out_path/ext_analysis/immuno.data.limma.asym0.prop90.Rda
limma_immuno_asym0_80: $out_path/ext_analysis/immuno.data.limma.asym0.prop80.Rda
limma_immuno_asym0_70: $out_path/ext_analysis/immuno.data.limma.asym0.prop70.Rda
limma_immuno_asym1_90: $out_path/ext_analysis/immuno.data.limma.asym1.prop90.Rda
limma_immuno_asym1_80: $out_path/ext_analysis/immuno.data.limma.asym1.prop80.Rda
limma_immuno_asym1_70: $out_path/ext_analysis/immuno.data.limma.asym1.prop70.Rda

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


# ALDEx3 gamma = 0 vs. 1: asymmetry around mean of normal distribution and proportion of genes NOT altered (thinning)
ALDEx3_0_immuno_asym0_90: $out_path/ext_analysis/immuno.data.aldex3_0.asym0.prop90.Rda
ALDEx3_0_immuno_asym0_80: $out_path/ext_analysis/immuno.data.aldex3_0.asym0.prop80.Rda
ALDEx3_0_immuno_asym0_70: $out_path/ext_analysis/immuno.data.aldex3_0.asym0.prop70.Rda
ALDEx3_1_immuno_asym0_90: $out_path/ext_analysis/immuno.data.aldex3_1.asym0.prop90.Rda
ALDEx3_1_immuno_asym0_80: $out_path/ext_analysis/immuno.data.aldex3_1.asym0.prop80.Rda
ALDEx3_1_immuno_asym0_70: $out_path/ext_analysis/immuno.data.aldex3_1.asym0.prop70.Rda

ALDEx3_0_immuno_asym1_90: $out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop90.Rda
ALDEx3_0_immuno_asym1_80: $out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop80.Rda
ALDEx3_0_immuno_asym1_70: $out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop70.Rda
ALDEx3_1_immuno_asym1_90: $out_path/ext_analysis/immuno.data.aldex3_1.asym1.prop90.Rda
ALDEx3_1_immuno_asym1_80: $out_path/ext_analysis/immuno.data.aldex3_1.asym1.prop80.Rda
ALDEx3_1_immuno_asym1_70: $out_path/ext_analysis/immuno.data.aldex3_1.asym1.prop70.Rda

ALDEx3_0_immuno_asym2_90: $out_path/ext_analysis/immuno.data.aldex3_0.asym2.prop90.Rda
ALDEx3_0_immuno_asym2_80: $out_path/ext_analysis/immuno.data.aldex3_0.asym2.prop80.Rda
ALDEx3_0_immuno_asym2_70: $out_path/ext_analysis/immuno.data.aldex3_0.asym2.prop70.Rda
ALDEx3_1_immuno_asym2_90: $out_path/ext_analysis/immuno.data.aldex3_1.asym2.prop90.Rda
ALDEx3_1_immuno_asym2_80: $out_path/ext_analysis/immuno.data.aldex3_1.asym2.prop80.Rda
ALDEx3_1_immuno_asym2_70: $out_path/ext_analysis/immuno.data.aldex3_1.asym2.prop70.Rda

ALDEx3_0_immuno_asym3_90: $out_path/ext_analysis/immuno.data.aldex3_0.asym3.prop90.Rda
ALDEx3_0_immuno_asym3_80: $out_path/ext_analysis/immuno.data.aldex3_0.asym3.prop80.Rda
ALDEx3_0_immuno_asym3_70: $out_path/ext_analysis/immuno.data.aldex3_0.asym3.prop70.Rda
ALDEx3_1_immuno_asym3_90: $out_path/ext_analysis/immuno.data.aldex3_1.asym3.prop90.Rda
ALDEx3_1_immuno_asym3_80: $out_path/ext_analysis/immuno.data.aldex3_1.asym3.prop80.Rda
ALDEx3_1_immuno_asym3_70: $out_path/ext_analysis/immuno.data.aldex3_1.asym3.prop70.Rda


### ALDEx3 gamma = 0 vs limma: asymmetry around mean of normal distribution and majority of genes ARE altered (thinning)
ALDEx3_0_immuno_asym1_49: $out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop49.Rda
limma_immuno_asym1_49: $out_path/ext_analysis/immuno.data.limma.asym1.prop49.Rda



# rules to generate the output files

### DESeq
$out_path/ext_analysis/brca.data.deseq.Rda : code/brca_deseq.R
	Rscript 'code/brca_deseq.R'

$out_path/ext_analysis/cdiff.data.deseq.Rda : code/cdiff_deseq.R 
	Rscript 'code/cdiff_deseq.R'

$out_path/ext_analysis/fre.data.deseq.Rda : code/fre_deseq.R 
	Rscript 'code/fre_deseq.R'

$out_path/ext_analysis/hiv.data.deseq.Rda : code/hiv_deseq.R 
	Rscript 'code/hiv_deseq.R'

$out_path/ext_analysis/immuno.data.deseq.Rda : code/immuno_deseq.R
	Rscript 'code/immuno_deseq.R'
	
$out_path/ext_analysis/kirc.data.deseq.Rda : code/kirc_deseq.R
	Rscript 'code/kirc_deseq.R'

$out_path/ext_analysis/lihc.data.deseq.Rda : code/lihc_deseq.R
	Rscript 'code/lihc_deseq.R'

$out_path/ext_analysis/luad.data.deseq.Rda : code/luad_deseq.R
	Rscript 'code/luad_deseq.R'

$out_path/ext_analysis/mar.data.deseq.Rda : code/mar_deseq.R 
	Rscript 'code/mar_deseq.R'
	
$out_path/ext_analysis/mts.data.deseq.Rda : code/mts_deseq.R
	Rscript 'code/mts_deseq.R'

$out_path/ext_analysis/oral.data.deseq.Rda : code/oral_deseq.R 
	Rscript 'code/oral_deseq.R'

$out_path/ext_analysis/prad.data.deseq.Rda : code/prad_deseq.R
	Rscript 'code/prad_deseq.R'

$out_path/ext_analysis/sccyto.data.deseq.Rda : code/sccyto_deseq.R 
	Rscript 'code/sccyto_deseq.R'

$out_path/ext_analysis/thca.data.deseq.Rda : code/thca_deseq.R
	Rscript 'code/thca_deseq.R'


### edgeR
$out_path/ext_analysis/brca.data.edger.Rda : code/brca_edgeR.R 
	Rscript 'code/brca_edgeR.R'

$out_path/ext_analysis/cdiff.data.edger.Rda : code/cdiff_edgeR.R
	Rscript 'code/cdiff_edgeR.R'
	
$out_path/ext_analysis/fre.data.edger.Rda : code/fre_edgeR.R
	Rscript 'code/fre_edgeR.R'

$out_path/ext_analysis/hiv.data.edger.Rda : code/hiv_edgeR.R
	Rscript 'code/hiv_edgeR.R'

$out_path/ext_analysis/immuno.data.edger.Rda : code/immuno_edgeR.R 
	Rscript 'code/immuno_edgeR.R'
	
$out_path/ext_analysis/kirc.data.edger.Rda : code/kirc_edgeR.R 
	Rscript 'code/kirc_edgeR.R'

$out_path/ext_analysis/lihc.data.edger.Rda : code/lihc_edgeR.R 
	Rscript 'code/lihc_edgeR.R'
	
$out_path/ext_analysis/luad.data.edger.Rda : code/luad_edgeR.R 
	Rscript 'code/luad_edgeR.R'
	
$out_path/ext_analysis/mar.data.edger.Rda : code/mar_edgeR.R
	Rscript 'code/mar_edgeR.R'

$out_path/ext_analysis/mts.data.edger.Rda : code/mts_edgeR.R
	Rscript 'code/mts_edgeR.R'

$out_path/ext_analysis/oral.data.edger.Rda : code/oral_edgeR.R
	Rscript 'code/oral_edgeR.R'

$out_path/ext_analysis/prad.data.edger.Rda : code/prad_edgeR.R 
	Rscript 'code/prad_edgeR.R'

$out_path/ext_analysis/sccyto.data.edger.Rda : code/sccyto_edgeR.R
	Rscript 'code/sccyto_edgeR.R'

$out_path/ext_analysis/thca.data.edger.Rda : code/thca_edgeR.R 
	Rscript 'code/thca_edgeR.R'


### limma 		
$out_path/ext_analysis/brca.data.limma.Rda : code/brca_limma.R 
	Rscript 'code/brca_limma.R'
	
$out_path/ext_analysis/cdiff.data.limma.Rda : code/cdiff_limma.R
	Rscript 'code/cdiff_limma.R'
	
$out_path/ext_analysis/fre.data.limma.Rda : code/fre_limma.R
	Rscript 'code/fre_limma.R'

$out_path/ext_analysis/hiv.data.limma.Rda : code/hiv_limma.R
	Rscript 'code/hiv_limma.R'

$out_path/ext_analysis/immuno.data.limma.Rda : code/immuno_limma.R 
	Rscript 'code/immuno_limma.R'

$out_path/ext_analysis/kirc.data.limma.Rda : code/kirc_limma.R 
	Rscript 'code/kirc_limma.R'
	
$out_path/ext_analysis/lihc.data.limma.Rda : code/lihc_limma.R 
	Rscript 'code/lihc_limma.R'

$out_path/ext_analysis/luad.data.limma.Rda : code/luad_limma.R 
	Rscript 'code/luad_limma.R'

$out_path/ext_analysis/mar.data.limma.Rda : code/mar_limma.R
	Rscript 'code/mar_limma.R'

$out_path/ext_analysis/mts.data.limma.Rda : code/mts_limma.R
	Rscript 'code/mts_limma.R'
	
$out_path/ext_analysis/oral.data.limma.Rda : code/oral_limma.R
	Rscript 'code/oral_limma.R'	

$out_path/ext_analysis/prad.data.limma.Rda : code/prad_limma.R 
	Rscript 'code/prad_limma.R'
	
$out_path/ext_analysis/sccyto.data.limma.Rda : code/sccyto_limma.R
	Rscript 'code/sccyto_limma.R'

$out_path/ext_analysis/thca.data.limma.Rda : code/thca_limma.R 
	Rscript 'code/thca_limma.R'


### aldex2_0
$out_path/ext_analysis/brca.data.aldex2_0.Rda : code/brca_aldex2_0.R 
	Rscript 'code/brca_aldex2_0.R'
	
$out_path/ext_analysis/cdiff.data.aldex2_0.Rda : code/cdiff_aldex2_0.R 
	Rscript 'code/cdiff_aldex2_0.R'
	
$out_path/ext_analysis/fre.data.aldex2_0.Rda : code/fre_aldex2_0.R 
	Rscript 'code/fre_aldex2_0.R'
	
$out_path/ext_analysis/hiv.data.aldex2_0.Rda : code/hiv_aldex2_0.R 
	Rscript 'code/hiv_aldex2_0.R'

$out_path/ext_analysis/immuno.data.aldex2_0.Rda : code/immuno_aldex2_0.R 
	Rscript 'code/immuno_aldex2_0.R'

$out_path/ext_analysis/kirc.data.aldex2_0.Rda : code/kirc_aldex2_0.R 
	Rscript 'code/kirc_aldex2_0.R'
	
$out_path/ext_analysis/lihc.data.aldex2_0.Rda : code/lihc_aldex2_0.R 
	Rscript 'code/lihc_aldex2_0.R'
	
$out_path/ext_analysis/luad.data.aldex2_0.Rda : code/luad_aldex2_0.R 
	Rscript 'code/luad_aldex2_0.R'
	
$out_path/ext_analysis/mar.data.aldex2_0.Rda : code/mar_aldex2_0.R 
	Rscript 'code/mar_aldex2_0.R'

$out_path/ext_analysis/mts.data.aldex2_0.Rda : code/mts_aldex2_0.R
	Rscript 'code/mts_aldex2_0.R'
	
$out_path/ext_analysis/oral.data.aldex2_0.Rda : code/oral_aldex2_0.R 
	Rscript 'code/oral_aldex2_0.R'

$out_path/ext_analysis/prad.data.aldex2_0.Rda : code/prad_aldex2_0.R 
	Rscript 'code/prad_aldex2_0.R'
	
$out_path/ext_analysis/sccyto.data.aldex2_0.Rda : code/sccyto_aldex2_0.R 
	Rscript 'code/sccyto_aldex2_0.R'

$out_path/ext_analysis/thca.data.aldex2_0.Rda : code/thca_aldex2_0.R 
	Rscript 'code/thca_aldex2_0.R'


### aldex2_1
$out_path/ext_analysis/brca.data.aldex2_1.Rda : code/brca_aldex2_1.R 
	Rscript 'code/brca_aldex2_1.R'
	
$out_path/ext_analysis/cdiff.data.aldex2_1.Rda : code/cdiff_aldex2_1.R 
	Rscript 'code/cdiff_aldex2_1.R'
	
$out_path/ext_analysis/fre.data.aldex2_1.Rda : code/fre_aldex2_1.R 
	Rscript 'code/fre_aldex2_1.R'
	
$out_path/ext_analysis/hiv.data.aldex2_1.Rda : code/hiv_aldex2_1.R 
	Rscript 'code/hiv_aldex2_1.R'

$out_path/ext_analysis/immuno.data.aldex2_1.Rda : code/immuno_aldex2_1.R
	Rscript 'code/immuno_aldex2_1.R'

$out_path/ext_analysis/kirc.data.aldex2_1.Rda : code/kirc_aldex2_1.R 
	Rscript 'code/kirc_aldex2_1.R'
	
$out_path/ext_analysis/lihc.data.aldex2_1.Rda : code/lihc_aldex2_1.R 
	Rscript 'code/lihc_aldex2_1.R'
	
$out_path/ext_analysis/luad.data.aldex2_1.Rda : code/luad_aldex2_1.R 
	Rscript 'code/luad_aldex2_1.R'
	
$out_path/ext_analysis/mar.data.aldex2_1.Rda : code/mar_aldex2_1.R 
	Rscript 'code/mar_aldex2_1.R'

$out_path/ext_analysis/mts.data.aldex2_1.Rda : code/mts_aldex2_1.R
	Rscript 'code/mts_aldex2_1.R'

$out_path/ext_analysis/oral.data.aldex2_1.Rda : code/oral_aldex2_1.R 
	Rscript 'code/oral_aldex2_1.R'

$out_path/ext_analysis/prad.data.aldex2_1.Rda : code/prad_aldex2_1.R 
	Rscript 'code/prad_aldex2_1.R'
	
$out_path/ext_analysis/sccyto.data.aldex2_1.Rda : code/sccyto_aldex2_1.R 
	Rscript 'code/sccyto_aldex2_1.R'

$out_path/ext_analysis/thca.data.aldex2_1.Rda : code/thca_aldex2_1.R 
	Rscript 'code/thca_aldex2_1.R'


### aldex2_2
$out_path/ext_analysis/brca.data.aldex2_2.Rda : code/brca_aldex2_2.R 
	Rscript 'code/brca_aldex2_2.R'
	
$out_path/ext_analysis/cdiff.data.aldex2_2.Rda : code/cdiff_aldex2_2.R 
	Rscript 'code/cdiff_aldex2_2.R'
	
$out_path/ext_analysis/fre.data.aldex2_2.Rda : code/fre_aldex2_2.R 
	Rscript 'code/fre_aldex2_2.R'
	
$out_path/ext_analysis/hiv.data.aldex2_2.Rda : code/hiv_aldex2_2.R 
	Rscript 'code/hiv_aldex2_2.R'

$out_path/ext_analysis/immuno.data.aldex2_2.Rda : code/immuno_aldex2_2.R 
	Rscript 'code/immuno_aldex2_2.R'

$out_path/ext_analysis/kirc.data.aldex2_2.Rda : code/kirc_aldex2_2.R 
	Rscript 'code/kirc_aldex2_2.R'
	
$out_path/ext_analysis/lihc.data.aldex2_2.Rda : code/lihc_aldex2_2.R 
	Rscript 'code/lihc_aldex2_2.R'
	
$out_path/ext_analysis/luad.data.aldex2_2.Rda : code/luad_aldex2_2.R 
	Rscript 'code/luad_aldex2_2.R'
	
$out_path/ext_analysis/mar.data.aldex2_2.Rda : code/mar_aldex2_2.R 
	Rscript 'code/mar_aldex2_2.R'

$out_path/ext_analysis/mts.data.aldex2_2.Rda : code/mts_aldex2_2.R
	Rscript 'code/mts_aldex2_2.R'

$out_path/ext_analysis/oral.data.aldex2_2.Rda : code/oral_aldex2_2.R 
	Rscript 'code/oral_aldex2_2.R'

$out_path/ext_analysis/prad.data.aldex2_2.Rda : code/prad_aldex2_2.R 
	Rscript 'code/prad_aldex2_2.R'
	
$out_path/ext_analysis/sccyto.data.aldex2_2.Rda : code/sccyto_aldex2_2.R 
	Rscript 'code/sccyto_aldex2_2.R'

$out_path/ext_analysis/thca.data.aldex2_2.Rda : code/thca_aldex2_2.R 
	Rscript 'code/thca_aldex2_2.R'


### aldex2_3
$out_path/ext_analysis/brca.data.aldex2_3.Rda : code/brca_aldex2_3.R 
	Rscript 'code/brca_aldex2_3.R'
	
$out_path/ext_analysis/cdiff.data.aldex2_3.Rda : code/cdiff_aldex2_3.R 
	Rscript 'code/cdiff_aldex2_3.R'
	
$out_path/ext_analysis/fre.data.aldex2_3.Rda : code/fre_aldex2_3.R 
	Rscript 'code/fre_aldex2_3.R'
	
$out_path/ext_analysis/hiv.data.aldex2_3.Rda : code/hiv_aldex2_3.R 
	Rscript 'code/hiv_aldex2_3.R'

$out_path/ext_analysis/immuno.data.aldex2_3.Rda : code/immuno_aldex2_3.R 
	Rscript 'code/immuno_aldex2_3.R'

$out_path/ext_analysis/kirc.data.aldex2_3.Rda : code/kirc_aldex2_3.R 
	Rscript 'code/kirc_aldex2_3.R'
	
$out_path/ext_analysis/lihc.data.aldex2_3.Rda : code/lihc_aldex2_3.R 
	Rscript 'code/lihc_aldex2_3.R'
	
$out_path/ext_analysis/luad.data.aldex2_3.Rda : code/luad_aldex2_3.R 
	Rscript 'code/luad_aldex2_3.R'
	
$out_path/ext_analysis/mar.data.aldex2_3.Rda : code/mar_aldex2_3.R 
	Rscript 'code/mar_aldex2_3.R'

$out_path/ext_analysis/mts.data.aldex2_3.Rda : code/mts_aldex2_3.R
	Rscript 'code/mts_aldex2_3.R'

$out_path/ext_analysis/oral.data.aldex2_3.Rda : code/oral_aldex2_3.R 
	Rscript 'code/oral_aldex2_3.R'

$out_path/ext_analysis/prad.data.aldex2_3.Rda : code/prad_aldex2_3.R 
	Rscript 'code/prad_aldex2_3.R'
	
$out_path/ext_analysis/sccyto.data.aldex2_3.Rda : code/sccyto_aldex2_3.R 
	Rscript 'code/sccyto_aldex2_3.R'

$out_path/ext_analysis/thca.data.aldex2_3.Rda : code/thca_aldex2_3.R 
	Rscript 'code/thca_aldex2_3.R'


### aldex2_4
$out_path/ext_analysis/brca.data.aldex2_4.Rda : code/brca_aldex2_4.R 
	Rscript 'code/brca_aldex2_4.R'
	
$out_path/ext_analysis/cdiff.data.aldex2_4.Rda : code/cdiff_aldex2_4.R 
	Rscript 'code/cdiff_aldex2_4.R'
	
$out_path/ext_analysis/fre.data.aldex2_4.Rda : code/fre_aldex2_4.R 
	Rscript 'code/fre_aldex2_4.R'
	
$out_path/ext_analysis/hiv.data.aldex2_4.Rda : code/hiv_aldex2_4.R 
	Rscript 'code/hiv_aldex2_4.R'

$out_path/ext_analysis/immuno.data.aldex2_4.Rda : code/immuno_aldex2_4.R 
	Rscript 'code/immuno_aldex2_4.R'

$out_path/ext_analysis/kirc.data.aldex2_4.Rda : code/kirc_aldex2_4.R 
	Rscript 'code/kirc_aldex2_4.R'
	
$out_path/ext_analysis/lihc.data.aldex2_4.Rda : code/lihc_aldex2_4.R 
	Rscript 'code/lihc_aldex2_4.R'
	
$out_path/ext_analysis/luad.data.aldex2_4.Rda : code/luad_aldex2_4.R 
	Rscript 'code/luad_aldex2_4.R'
	
$out_path/ext_analysis/mar.data.aldex2_4.Rda : code/mar_aldex2_4.R 
	Rscript 'code/mar_aldex2_4.R'

$out_path/ext_analysis/mts.data.aldex2_4.Rda : code/mts_aldex2_4.R
	Rscript 'code/mts_aldex2_4.R'

$out_path/ext_analysis/oral.data.aldex2_4.Rda : code/oral_aldex2_4.R 
	Rscript 'code/oral_aldex2_4.R'

$out_path/ext_analysis/prad.data.aldex2_4.Rda : code/prad_aldex2_4.R 
	Rscript 'code/prad_aldex2_4.R'
	
$out_path/ext_analysis/sccyto.data.aldex2_4.Rda : code/sccyto_aldex2_4.R 
	Rscript 'code/sccyto_aldex2_4.R'

$out_path/ext_analysis/thca.data.aldex2_4.Rda : code/thca_aldex2_4.R 
	Rscript 'code/thca_aldex2_4.R'


### aldex2_5
$out_path/ext_analysis/brca.data.aldex2_5.Rda : code/brca_aldex2_5.R 
	Rscript 'code/brca_aldex2_5.R'
	
$out_path/ext_analysis/cdiff.data.aldex2_5.Rda : code/cdiff_aldex2_5.R 
	Rscript 'code/cdiff_aldex2_5.R'
	
$out_path/ext_analysis/fre.data.aldex2_5.Rda : code/fre_aldex2_5.R 
	Rscript 'code/fre_aldex2_5.R'

$out_path/ext_analysis/hiv.data.aldex2_5.Rda : code/hiv_aldex2_5.R 
	Rscript 'code/hiv_aldex2_5.R'

$out_path/ext_analysis/immuno.data.aldex2_5.Rda : code/immuno_aldex2_5.R 
	Rscript 'code/immuno_aldex2_5.R'

$out_path/ext_analysis/kirc.data.aldex2_5.Rda : code/kirc_aldex2_5.R 
	Rscript 'code/kirc_aldex2_5.R'
	
$out_path/ext_analysis/lihc.data.aldex2_5.Rda : code/lihc_aldex2_5.R 
	Rscript 'code/lihc_aldex2_5.R'
	
$out_path/ext_analysis/luad.data.aldex2_5.Rda : code/luad_aldex2_5.R 
	Rscript 'code/luad_aldex2_5.R'
	
$out_path/ext_analysis/mar.data.aldex2_5.Rda : code/mar_aldex2_5.R 
	Rscript 'code/mar_aldex2_5.R'

$out_path/ext_analysis/mts.data.aldex2_5.Rda : code/mts_aldex2_5.R
	Rscript 'code/mts_aldex2_5.R'

$out_path/ext_analysis/oral.data.aldex2_5.Rda : code/oral_aldex2_5.R 
	Rscript 'code/oral_aldex2_5.R'

$out_path/ext_analysis/prad.data.aldex2_5.Rda : code/prad_aldex2_5.R 
	Rscript 'code/prad_aldex2_5.R'
	
$out_path/ext_analysis/sccyto.data.aldex2_5.Rda : code/sccyto_aldex2_5.R 
	Rscript 'code/sccyto_aldex2_5.R'

$out_path/ext_analysis/thca.data.aldex2_5.Rda : code/thca_aldex2_5.R 
	Rscript 'code/thca_aldex2_5.R'


### aldex3_0
$out_path/ext_analysis/brca.data.aldex3_0.Rda : code/brca_aldex3_0.R 
	Rscript 'code/brca_aldex3_0.R'
	
$out_path/ext_analysis/cdiff.data.aldex3_0.Rda : code/cdiff_aldex3_0.R 
	Rscript 'code/cdiff_aldex3_0.R'
	
$out_path/ext_analysis/fre.data.aldex3_0.Rda : code/fre_aldex3_0.R 
	Rscript 'code/fre_aldex3_0.R'
	
$out_path/ext_analysis/hiv.data.aldex3_0.Rda : code/hiv_aldex3_0.R 
	Rscript 'code/hiv_aldex3_0.R'

$out_path/ext_analysis/immuno.data.aldex3_0.Rda : code/immuno_aldex3_0.R 
	Rscript 'code/immuno_aldex3_0.R'

$out_path/ext_analysis/kirc.data.aldex3_0.Rda : code/kirc_aldex3_0.R 
	Rscript 'code/kirc_aldex3_0.R'
	
$out_path/ext_analysis/lihc.data.aldex3_0.Rda : code/lihc_aldex3_0.R 
	Rscript 'code/lihc_aldex3_0.R'
	
$out_path/ext_analysis/luad.data.aldex3_0.Rda : code/luad_aldex3_0.R 
	Rscript 'code/luad_aldex3_0.R'
	
$out_path/ext_analysis/mar.data.aldex3_0.Rda : code/mar_aldex3_0.R 
	Rscript 'code/mar_aldex3_0.R'

$out_path/ext_analysis/mts.data.aldex3_0.Rda : code/mts_aldex3_0.R
	Rscript 'code/mts_aldex3_0.R'
	
$out_path/ext_analysis/oral.data.aldex3_0.Rda : code/oral_aldex3_0.R 
	Rscript 'code/oral_aldex3_0.R'

$out_path/ext_analysis/prad.data.aldex3_0.Rda : code/prad_aldex3_0.R 
	Rscript 'code/prad_aldex3_0.R'
	
$out_path/ext_analysis/sccyto.data.aldex3_0.Rda : code/sccyto_aldex3_0.R 
	Rscript 'code/sccyto_aldex3_0.R'

$out_path/ext_analysis/thca.data.aldex3_0.Rda : code/thca_aldex3_0.R 
	Rscript 'code/thca_aldex3_0.R'


### aldex3_1
$out_path/ext_analysis/brca.data.aldex3_1.Rda : code/brca_aldex3_1.R 
	Rscript 'code/brca_aldex3_1.R'
	
$out_path/ext_analysis/cdiff.data.aldex3_1.Rda : code/cdiff_aldex3_1.R 
	Rscript 'code/cdiff_aldex3_1.R'
	
$out_path/ext_analysis/fre.data.aldex3_1.Rda : code/fre_aldex3_1.R 
	Rscript 'code/fre_aldex3_1.R'
	
$out_path/ext_analysis/hiv.data.aldex3_1.Rda : code/hiv_aldex3_1.R 
	Rscript 'code/hiv_aldex3_1.R'

$out_path/ext_analysis/immuno.data.aldex3_1.Rda : code/immuno_aldex3_1.R
	Rscript 'code/immuno_aldex3_1.R'

$out_path/ext_analysis/kirc.data.aldex3_1.Rda : code/kirc_aldex3_1.R 
	Rscript 'code/kirc_aldex3_1.R'
	
$out_path/ext_analysis/lihc.data.aldex3_1.Rda : code/lihc_aldex3_1.R 
	Rscript 'code/lihc_aldex3_1.R'
	
$out_path/ext_analysis/luad.data.aldex3_1.Rda : code/luad_aldex3_1.R 
	Rscript 'code/luad_aldex3_1.R'
	
$out_path/ext_analysis/mar.data.aldex3_1.Rda : code/mar_aldex3_1.R 
	Rscript 'code/mar_aldex3_1.R'

$out_path/ext_analysis/mts.data.aldex3_1.Rda : code/mts_aldex3_1.R
	Rscript 'code/mts_aldex3_1.R'

$out_path/ext_analysis/oral.data.aldex3_1.Rda : code/oral_aldex3_1.R 
	Rscript 'code/oral_aldex3_1.R'

$out_path/ext_analysis/prad.data.aldex3_1.Rda : code/prad_aldex3_1.R 
	Rscript 'code/prad_aldex3_1.R'
	
$out_path/ext_analysis/sccyto.data.aldex3_1.Rda : code/sccyto_aldex3_1.R 
	Rscript 'code/sccyto_aldex3_1.R'

$out_path/ext_analysis/thca.data.aldex3_1.Rda : code/thca_aldex3_1.R 
	Rscript 'code/thca_aldex3_1.R'


### aldex3_2
$out_path/ext_analysis/brca.data.aldex3_2.Rda : code/brca_aldex3_2.R 
	Rscript 'code/brca_aldex3_2.R'
	
$out_path/ext_analysis/cdiff.data.aldex3_2.Rda : code/cdiff_aldex3_2.R 
	Rscript 'code/cdiff_aldex3_2.R'
	
$out_path/ext_analysis/fre.data.aldex3_2.Rda : code/fre_aldex3_2.R 
	Rscript 'code/fre_aldex3_2.R'
	
$out_path/ext_analysis/hiv.data.aldex3_2.Rda : code/hiv_aldex3_2.R 
	Rscript 'code/hiv_aldex3_2.R'

$out_path/ext_analysis/immuno.data.aldex3_2.Rda : code/immuno_aldex3_2.R 
	Rscript 'code/immuno_aldex3_2.R'

$out_path/ext_analysis/kirc.data.aldex3_2.Rda : code/kirc_aldex3_2.R 
	Rscript 'code/kirc_aldex3_2.R'
	
$out_path/ext_analysis/lihc.data.aldex3_2.Rda : code/lihc_aldex3_2.R 
	Rscript 'code/lihc_aldex3_2.R'
	
$out_path/ext_analysis/luad.data.aldex3_2.Rda : code/luad_aldex3_2.R 
	Rscript 'code/luad_aldex3_2.R'
	
$out_path/ext_analysis/mar.data.aldex3_2.Rda : code/mar_aldex3_2.R 
	Rscript 'code/mar_aldex3_2.R'

$out_path/ext_analysis/mts.data.aldex3_2.Rda : code/mts_aldex3_2.R
	Rscript 'code/mts_aldex3_2.R'

$out_path/ext_analysis/oral.data.aldex3_2.Rda : code/oral_aldex3_2.R 
	Rscript 'code/oral_aldex3_2.R'

$out_path/ext_analysis/prad.data.aldex3_2.Rda : code/prad_aldex3_2.R 
	Rscript 'code/prad_aldex3_2.R'
	
$out_path/ext_analysis/sccyto.data.aldex3_2.Rda : code/sccyto_aldex3_2.R 
	Rscript 'code/sccyto_aldex3_2.R'

$out_path/ext_analysis/thca.data.aldex3_2.Rda : code/thca_aldex3_2.R 
	Rscript 'code/thca_aldex3_2.R'


### aldex3_3
$out_path/ext_analysis/brca.data.aldex3_3.Rda : code/brca_aldex3_3.R 
	Rscript 'code/brca_aldex3_3.R'
	
$out_path/ext_analysis/cdiff.data.aldex3_3.Rda : code/cdiff_aldex3_3.R 
	Rscript 'code/cdiff_aldex3_3.R'
	
$out_path/ext_analysis/fre.data.aldex3_3.Rda : code/fre_aldex3_3.R 
	Rscript 'code/fre_aldex3_3.R'
	
$out_path/ext_analysis/hiv.data.aldex3_3.Rda : code/hiv_aldex3_3.R 
	Rscript 'code/hiv_aldex3_3.R'

$out_path/ext_analysis/immuno.data.aldex3_3.Rda : code/immuno_aldex3_3.R 
	Rscript 'code/immuno_aldex3_3.R'

$out_path/ext_analysis/kirc.data.aldex3_3.Rda : code/kirc_aldex3_3.R 
	Rscript 'code/kirc_aldex3_3.R'
	
$out_path/ext_analysis/lihc.data.aldex3_3.Rda : code/lihc_aldex3_3.R 
	Rscript 'code/lihc_aldex3_3.R'
	
$out_path/ext_analysis/luad.data.aldex3_3.Rda : code/luad_aldex3_3.R 
	Rscript 'code/luad_aldex3_3.R'
	
$out_path/ext_analysis/mar.data.aldex3_3.Rda : code/mar_aldex3_3.R 
	Rscript 'code/mar_aldex3_3.R'

$out_path/ext_analysis/mts.data.aldex3_3.Rda : code/mts_aldex3_3.R
	Rscript 'code/mts_aldex3_3.R'

$out_path/ext_analysis/oral.data.aldex3_3.Rda : code/oral_aldex3_3.R 
	Rscript 'code/oral_aldex3_3.R'

$out_path/ext_analysis/prad.data.aldex3_3.Rda : code/prad_aldex3_3.R 
	Rscript 'code/prad_aldex3_3.R'
	
$out_path/ext_analysis/sccyto.data.aldex3_3.Rda : code/sccyto_aldex3_3.R 
	Rscript 'code/sccyto_aldex3_3.R'

$out_path/ext_analysis/thca.data.aldex3_3.Rda : code/thca_aldex3_3.R 
	Rscript 'code/thca_aldex3_3.R'


### aldex3_4
$out_path/ext_analysis/brca.data.aldex3_4.Rda : code/brca_aldex3_4.R 
	Rscript 'code/brca_aldex3_4.R'
	
$out_path/ext_analysis/cdiff.data.aldex3_4.Rda : code/cdiff_aldex3_4.R 
	Rscript 'code/cdiff_aldex3_4.R'
	
$out_path/ext_analysis/fre.data.aldex3_4.Rda : code/fre_aldex3_4.R 
	Rscript 'code/fre_aldex3_4.R'
	
$out_path/ext_analysis/hiv.data.aldex3_4.Rda : code/hiv_aldex3_4.R 
	Rscript 'code/hiv_aldex3_4.R'

$out_path/ext_analysis/immuno.data.aldex3_4.Rda : code/immuno_aldex3_4.R 
	Rscript 'code/immuno_aldex3_4.R'

$out_path/ext_analysis/kirc.data.aldex3_4.Rda : code/kirc_aldex3_4.R 
	Rscript 'code/kirc_aldex3_4.R'
	
$out_path/ext_analysis/lihc.data.aldex3_4.Rda : code/lihc_aldex3_4.R 
	Rscript 'code/lihc_aldex3_4.R'
	
$out_path/ext_analysis/luad.data.aldex3_4.Rda : code/luad_aldex3_4.R 
	Rscript 'code/luad_aldex3_4.R'
	
$out_path/ext_analysis/mar.data.aldex3_4.Rda : code/mar_aldex3_4.R 
	Rscript 'code/mar_aldex3_4.R'

$out_path/ext_analysis/mts.data.aldex3_4.Rda : code/mts_aldex3_4.R
	Rscript 'code/mts_aldex3_4.R'

$out_path/ext_analysis/oral.data.aldex3_4.Rda : code/oral_aldex3_4.R 
	Rscript 'code/oral_aldex3_4.R'

$out_path/ext_analysis/prad.data.aldex3_4.Rda : code/prad_aldex3_4.R 
	Rscript 'code/prad_aldex3_4.R'
	
$out_path/ext_analysis/sccyto.data.aldex3_4.Rda : code/sccyto_aldex3_4.R 
	Rscript 'code/sccyto_aldex3_4.R'

$out_path/ext_analysis/thca.data.aldex3_4.Rda : code/thca_aldex3_4.R 
	Rscript 'code/thca_aldex3_4.R'


### aldex3_5
$out_path/ext_analysis/brca.data.aldex3_5.Rda : code/brca_aldex3_5.R 
	Rscript 'code/brca_aldex3_5.R'
	
$out_path/ext_analysis/cdiff.data.aldex3_5.Rda : code/cdiff_aldex3_5.R 
	Rscript 'code/cdiff_aldex3_5.R'
	
$out_path/ext_analysis/fre.data.aldex3_5.Rda : code/fre_aldex3_5.R 
	Rscript 'code/fre_aldex3_5.R'

$out_path/ext_analysis/hiv.data.aldex3_5.Rda : code/hiv_aldex3_5.R 
	Rscript 'code/hiv_aldex3_5.R'

$out_path/ext_analysis/immuno.data.aldex3_5.Rda : code/immuno_aldex3_5.R 
	Rscript 'code/immuno_aldex3_5.R'

$out_path/ext_analysis/kirc.data.aldex3_5.Rda : code/kirc_aldex3_5.R 
	Rscript 'code/kirc_aldex3_5.R'
	
$out_path/ext_analysis/lihc.data.aldex3_5.Rda : code/lihc_aldex3_5.R 
	Rscript 'code/lihc_aldex3_5.R'
	
$out_path/ext_analysis/luad.data.aldex3_5.Rda : code/luad_aldex3_5.R 
	Rscript 'code/luad_aldex3_5.R'
	
$out_path/ext_analysis/mar.data.aldex3_5.Rda : code/mar_aldex3_5.R 
	Rscript 'code/mar_aldex3_5.R'

$out_path/ext_analysis/mts.data.aldex3_5.Rda : code/mts_aldex3_5.R
	Rscript 'code/mts_aldex3_5.R'

$out_path/ext_analysis/oral.data.aldex3_5.Rda : code/oral_aldex3_5.R 
	Rscript 'code/oral_aldex3_5.R'

$out_path/ext_analysis/prad.data.aldex3_5.Rda : code/prad_aldex3_5.R 
	Rscript 'code/prad_aldex3_5.R'
	
$out_path/ext_analysis/sccyto.data.aldex3_5.Rda : code/sccyto_aldex3_5.R 
	Rscript 'code/sccyto_aldex3_5.R'

$out_path/ext_analysis/thca.data.aldex3_5.Rda : code/thca_aldex3_5.R 
	Rscript 'code/thca_aldex3_5.R'


### limma: asymmetry around mean of normal distribution (thinning)
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


### limma: asymmetry around mean of normal distribution and proportion of genes NOT altered (thinning)
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


### ALDEx3 gamma = 0 vs. 1: asymmetry around mean of normal distribution and proportion of genes NOT altered (thinning)
$out_path/ext_analysis/immuno.data.aldex3_0.asym0.prop90.Rda : code/immuno_aldex3_0_asym0_90.R
	Rscript 'code/immuno_aldex3_0_asym0_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym0.prop80.Rda : code/immuno_aldex3_0_asym0_80.R
	Rscript 'code/immuno_aldex3_0_asym0_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym0.prop70.Rda : code/immuno_aldex3_0_asym0_70.R
	Rscript 'code/immuno_aldex3_0_asym0_70.R'

$out_path/ext_analysis/immuno.data.aldex3_1.asym0.prop90.Rda  : code/immuno_aldex3_1_asym0_90.R
	Rscript 'code/immuno_aldex3_1_asym0_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym0.prop80.Rda : code/immuno_aldex3_1_asym0_80.R
	Rscript 'code/immuno_aldex3_1_asym0_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym0.prop70.Rda : code/immuno_aldex3_1_asym0_70.R
	Rscript 'code/immuno_aldex3_1_asym0_70.R'

$out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop90.Rda  : code/immuno_aldex3_0_asym1_90.R
	Rscript 'code/immuno_aldex3_0_asym1_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop80.Rda  : code/immuno_aldex3_0_asym1_80.R
	Rscript 'code/immuno_aldex3_0_asym1_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop70.Rda  : code/immuno_aldex3_0_asym1_70.R
	Rscript 'code/immuno_aldex3_0_asym1_70.R'

$out_path/ext_analysis/immuno.data.aldex3_1.asym1.prop90.Rda  : code/immuno_aldex3_1_asym1_90.R
	Rscript 'code/immuno_aldex3_1_asym1_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym1.prop80.Rda  : code/immuno_aldex3_1_asym1_80.R
	Rscript 'code/immuno_aldex3_1_asym1_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym1.prop70.Rda  : code/immuno_aldex3_1_asym1_70.R
	Rscript 'code/immuno_aldex3_1_asym1_70.R'

$out_path/ext_analysis/immuno.data.aldex3_0.asym2.prop90.Rda : code/immuno_aldex3_0_asym2_90.R
	Rscript 'code/immuno_aldex3_0_asym2_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym2.prop80.Rda : code/immuno_aldex3_0_asym2_80.R
	Rscript 'code/immuno_aldex3_0_asym2_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym2.prop70.Rda : code/immuno_aldex3_0_asym2_70.R
	Rscript 'code/immuno_aldex3_0_asym2_70.R'

$out_path/ext_analysis/immuno.data.aldex3_1.asym2.prop90.Rda : code/immuno_aldex3_1_asym2_90.R
	Rscript 'code/immuno_aldex3_1_asym2_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym2.prop80.Rda : code/immuno_aldex3_1_asym2_80.R
	Rscript 'code/immuno_aldex3_1_asym2_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym2.prop70.Rda : code/immuno_aldex3_1_asym2_70.R
	Rscript 'code/immuno_aldex3_1_asym2_70.R'

$out_path/ext_analysis/immuno.data.aldex3_0.asym3.prop90.Rda : code/immuno_aldex3_0_asym3_90.R
	Rscript 'code/immuno_aldex3_0_asym3_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym3.prop80.Rda : code/immuno_aldex3_0_asym3_80.R
	Rscript 'code/immuno_aldex3_0_asym3_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_0.asym3.prop70.Rda : code/immuno_aldex3_0_asym3_70.R
	Rscript 'code/immuno_aldex3_0_asym3_70.R'

$out_path/ext_analysis/immuno.data.aldex3_1.asym3.prop90.Rda : code/immuno_aldex3_1_asym3_90.R
	Rscript 'code/immuno_aldex3_1_asym3_90.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym3.prop80.Rda : code/immuno_aldex3_1_asym3_80.R
	Rscript 'code/immuno_aldex3_1_asym3_80.R'
	
$out_path/ext_analysis/immuno.data.aldex3_1.asym3.prop70.Rda : code/immuno_aldex3_1_asym3_70.R
	Rscript 'code/immuno_aldex3_1_asym3_70.R'


### ALDEx3 gamma = 0 vs limma: asymmetry around mean of normal distribution and majority of genes ARE altered (thinning)
$out_path/ext_analysis/immuno.data.aldex3_0.asym1.prop49.Rda : code/immuno_aldex3_0_asym1_49.R
	Rscript 'code/immuno_aldex3_0_asym1_49.R'

$out_path/ext_analysis/immuno.data.limma.asym1.prop49.Rda : code/immuno_limma_asym1_49.R
	Rscript 'code/immuno_limma_asym1_49.R'



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

	
	
