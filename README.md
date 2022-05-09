# scRNA analysis of hematopoiesis in aging and disease

This repository contains the scripts used for the analysis shown in **Uncovering perturbations in human hematopoiesis associated with healthy aging and myeloid malignancies at single cell resolution**.

Folders contain the following scripts: 

**01_integration:** integration of samples and unsupervised clustering using Seurat
  - 01_explore_individual_sample: creation of Seurat objects from individual raw count matrices created with CellRanger, QC filtering and exploratory plots.
  - 02a_integrate_samples_young: integration of 5 young samples, unsupervised clustering and manual annotation. 
  - 02b_integrate_samples_senior: integration of 3 elderly samples.
  - 02c_integrate_samples_MDS: integration of 4 MDS samples.
  - 03_integrate_samples_different_condition: integration of young and elderly to create a shared UMAP
  - 04_proportion_test: test for differences in cell type proportion

**02_glmnet_classification:** scripts for the cell type classification method based on GLMnet
  - 01_binary_models: build classification models for individual cell types
  - 02_final_classification: assign final cell type labels by comparing the results from the binary models and choosing the one with higher scores.

**03_differential_expression:** differential expression analysis between cell types and conditions and subsequent GSEA
  - 01_differential_expression: script to perform differential expression
  - 02_GSEA_young_elderly: GSEA for differential expression between young and elderly and plot results
  - 03_GSEA_young_elderly_mds: GSEA for differential expression between MDS and both young and elderly and plot results

**04_trajectory_analysis:** scripts to perform trajectory inference with Stream and Palantir and downstream analysis
  - **01_stream:** scripts for Stream
      - GenerateData: prepare data to run Stream
      - Stream: run Stream on young samples and project elderly samples on the resulting trajectory
  - **02_palantir:** scripts for Palantir  
