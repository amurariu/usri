# Code and study data for Dos Santos et al. (2025)

Last updated: 28th October 2025

This repository contains input data and code for running all analyses and producing all figures in the article, **"Are any of these results even real?!": one simple rule to control the false discovery rate when analysing high-throughput sequencing data**, published in XYZ.

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

