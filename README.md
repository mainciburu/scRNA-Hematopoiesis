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

