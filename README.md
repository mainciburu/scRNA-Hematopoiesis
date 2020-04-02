# scRNA MDS project

01_explore_individual_sample => filtering, normalization, scaling, PCA, UMAP and clustering of a single sample 10X matrix
02_integrate_samples_same_condition.r => integration, scaling, PCA, UMAP and clustering of multiple preprocessed (filtered, normalized) seurat objects sharing the same condition
03_integrate_samples_different_condition.r => integration, scaling, PCA and UMAP of multiple preprocessed (filtered, normalized) and possibly already integrated seurat objects from different conditions
04_glmnet_classification/ => scripts for building glmnet binary models and predicting cell types in new data
	celltypes.txt => list of cell types to build models on
	SendArray.sh
	SlurmBatch_binary_models_R.sbs
	binary_models.r
	SlurmBatch_final_classification_R_final.sbs
	final_classification.r