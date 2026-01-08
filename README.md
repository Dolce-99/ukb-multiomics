# Multi-omics integration predicts 17 disease incidences in the UK Biobank

This repo contains all scripts for analysis pipeline in paper "Multi-omics integration predicts 17 disease incidences in the UK Biobank". 

[Figure 1.pdf](https://github.com/user-attachments/files/24510607/Figure.1.pdf)

Leveraging data from 23,776 UK Biobank participants, we integrated 159 NMR-based metabolites and 2,923 Olink affinity-based proteins to assess their incremental value over 3 traditional clinical predictor sets. 
We first fit Cox Proportional Hazard models for each baseline predictor sets to adjust out the covariate effect and obtained martingale residuals. Then we did omics integration to predict these residuals.
Beyond performance benchmarking (evalutaed by C-index), we further assessed the top contributing omics features, conducted KEGG pathway analysis, and validated potential causal relationships using Mendelian Randomization.
