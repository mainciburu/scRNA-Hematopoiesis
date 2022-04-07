import numpy as np
import pandas as pd     
import pickle
import matplotlib
import matplotlib.pyplot as plt
import scipy
import statsmodels
import palantir
import scanpy as sc
import seaborn as sns
from collections import Counter

######### Input data ###################

# loom/h5ad object, to extract meta data
h5ad_file="/home/mainciburu/scRNA/palantir/processed/senior_int.h5ad"
loom_file="/home/mainciburu/scRNA/palantir/data/senior_int.loom"
# Integrated expression matrix (all common genes)
int_df_file="/home/mainciburu/scRNA/palantir/data/senior_int_mat.csv"
# Normalized expression matrix (all common genes)
norm_df_full_file="/home/mainciburu/scRNA/palantir/data/senior_norm_mat_full.csv"
# UMAP coordinates
umap_seurat_file="/home/mainciburu/scRNA/palantir/data/senior_umap.txt"
# Names for final cells (from knn)
final_cells_file="/home/mainciburu/scRNA/palantir/data/final_cells_senior.csv"

plot_name = "/home/mainciburu/scRNA/palantir/pics/senior"
res_path="/home/mainciburu/scRNA/palantir/results/senior/"

#######################################

# Load loom object or h5ad processed object
adata = sc.read_loom(loom_file)
adata.write(filename=h5ad_file)
adata = sc.read(filename=h5ad_file)

# Integrated or RNA matrix from csv
int_df=pd.read_csv(int_df_file, index_col=0)
# transpose
int_df=int_df.T

# Project results on seurat umap
umap_seurat=pd.read_csv(umap_seurat_file, sep="\t", header=None, index_col = 0)
umap_seurat.columns = ['x', 'y']

# Preprocess data for Palantir
np.random.seed(123)
# PCA
pca_projections, _ = palantir.utils.run_pca(int_df, use_hvg=False)
# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=30)
# Multiscale data
ms_data = palantir.utils.determine_multiscale_space(dm_res)


### Select the cells with higher hsc probability
df = adata.obs['prediction']
df = pd.DataFrame(df)
df['prediction_prob'] = adata.obs['prediction_prob']
df = df.loc[df['prediction'] == "HSC"]
df = df.sort_values('prediction_prob', ascending = False)
start_cells = df.index[0]

# Import final cells from knn
final_cells = pd.read_csv(final_cells_file, header=0)
final_cells = list(final_cells['x'])

# Trajectories with ending cells and highest HSC
ms_data.index = adata.obs.index
pr_res = palantir.core.run_palantir(ms_data, start_cells[0], terminal_states=final_cells, n_jobs = 1)

# Rename trajectories
i=pr_res.branch_probs.columns
df = adata.obs['prediction']
df = df.loc[i]
names = list(df)
pr_res.branch_probs.columns = names

# Plot results
palantir.plot.highlight_cells_on_tsne(umap_seurat, final_cells)
plt.savefig(plot_name + "_umap_final_cells_knn.pdf")

palantir.plot.plot_palantir_results(pr_res, umap_seurat)
plt.savefig(plot_name + "_umap_branches_knn.pdf")

### Store results ###

# pseudotime, branch probability, differentiation potential
pr_res.branch_probs.to_csv(res_path + "branch_probs.csv")
pr_res.pseudotime.to_csv(res_path + "pseudotime.csv", )
pr_res.entropy.to_csv(res_path + "diff_potential.csv")

# Save pr_res object
file_pr = open(res_path + 'pr_res.obj', 'wb') 
pickle.dump(pr_res, file_pr)

# Imputed data for gene expression trends
# !!Use original normalized expression values (not integrated)
norm_df_full = pd.read_csv(norm_df_full_file, index_col = 0)
norm_df_full = pd.DataFrame.transpose(norm_df_full)
# Impute data
imp_df = palantir.utils.run_magic_imputation(norm_df_full, dm_res)
imp_df.to_hdf(res_path + "imp_df.hdf5", mode="w", key="imp_df")
