# Code and study data for Dos Santos et al. (2026)

Last updated: 29th May 2026

This repository contains input data and code for running all analyses and producing all figures in the article, **Not every gene is special: modelling scale controls the false discovery rate when analysing high-throughput sequencing data**, submitted to PLoS Computational Biology.

### Contents:
- **analysis:**
  - Directory containing the '.Rda' objects for:
    - Confusion matrices (i.e. the output of 'code/get_confusion.R').
    - Summary stats (contains minimal difference between groups needed for *P* <0.05, for ALDEx2 and ALDEx3).
   
- **code:**
  - Directory containing all R scripts for producing the figures in the paper and running all tools on all datasets w/ binomial thinning, over 100 iterations of the data.

- **data:**
  - Directory containing input & metadata for all datasets, as well as metadata and any functional lookup data.

- **figures:**
  - Directory containing .png files for all figures produced by 'code/fig_xyz' scripts.
  - Also contains sub-directory 'SVGs' containing the SVG files of the final figures used in the manuscript, and 'figureEdits.txt', describing exactly how all R-produced figures have been edited in Inkscape, if applicable.

- **paper:**
  - Directory containing all files relating to the initial submission of this study to PLOS Computational Biology.
  - Directory containing all files relating to the revised submission of this study to PLOS Computational Biology, including the supplementary material and latexdiff files showing marked up changes from the original submission.
  - Directory containing the final figures used in the manuscript main text and supplementary material.
