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
      - 01_seurat_to_loom: prepare data to run Palantir
      - 02_palantir_young: run Palantir on young samples
      - 03_knn_final_cells: find knn cells in elderly samples to use as final states in Palantir
      - 04_palantir_elderly: run Palantir on elderly samples
      - 05_palantir_mds: run Palantir on elderly samples
      - 06_palantir_stats: test for differences in Palantir results between young and elderly
      - 07_compute_gene_trends_script: run Palantir gene trends 
      - 08_read_trends_and_cluster: cluster gene trends
      - 09_monocyte_analysis: downstream analysis for the comparison of monocytes branch in young and elderly
      - 10_erythroid_analysis: downstream analysis for the comparison of erythroid branch in young, elderly and MDS

**05_GRN:** gene regulatory networks analysis
  - 01_GenerateData: prepare data for scenic
  - 02_pyscenic: run python implementation of scenic
  - 03_RSS: calculate regulon specificity score per cell type
  - 04_downstream_analysis: create regulons heatmap based on scenic results and perform term over-representation analysis
  - 05_CytoscapeVisualization: format scenic results for visualization in cytoscape

**figures:** additional scripts to reproduce paper figures.
