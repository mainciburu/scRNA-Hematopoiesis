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
h5ad_file="/home/mainciburu/scRNA/palantir/processed/mds5.h5ad"
loom_file="/home/mainciburu/scRNA/palantir/data/mds5_norm.loom"

# Normalized expression matrix (all common genes)
norm_df_full_file="/home/mainciburu/scRNA/palantir/data/mds5_norm_mat_full.csv"
# UMAP coordinates
umap_seurat_file="/home/mainciburu/scRNA/palantir/data/mds5_umap.txt"

plot_name = "/home/mainciburu/scRNA/palantir/pics/mds5"
res_path="/home/mainciburu/scRNA/palantir/results/mds5/"

#######################################

# Load loom object or h5ad processed object

adata = sc.read_loom(loom_file)
adata.write(filename=h5ad_file)
adata = sc.read(filename=h5ad_file)

# RNA matrix from csv
norm_df=pd.read_csv(norm_df_full_file, index_col=0)
# transpose
norm_df=norm_df.T

# Load seurat umap
umap_seurat=pd.read_csv(umap_seurat_file, sep="\t", header=None, index_col = 0)
umap_seurat.columns = ['x', 'y']

# Preprocess data for Palantir
np.random.seed(123)
# PCA
pca_projections, _ = palantir.utils.run_pca(norm_df, use_hvg=False)
# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=30)
# plot dm components on umap
dm_res["EigenVectors"].index = adata.obs.index
palantir.plot.plot_diffusion_components(umap_seurat, dm_res)
plt.savefig(plot_name + "_umap_dm_components.pdf")

# Multiscale data
ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=6)

### Select the cells with higher hsc probability
df = adata.obs['prediction']
df = pd.DataFrame(df)
df['prediction_prob'] = adata.obs['prediction_prob']
df = df.loc[df['prediction'] == "HSC"]
df = df.sort_values('prediction_prob', ascending = False)
start_cells = df.index[0:10]
palantir.plot.highlight_cells_on_tsne(umap_seurat, start_cells)
plt.savefig(plot_name + "_umap_start_cells.pdf")

### Compute trajectories for each start cell
end_cells = []
ms_data.index = adata.obs.index
for cell in start_cells:
    print("Starting " + cell + " -------")
    pr_res = palantir.core.run_palantir(ms_data, cell, n_jobs = 1)
    end_cells.append(pr_res.branch_probs.columns)

# Count results
end_cells2=[]
for i in end_cells:
    for j in i:
        end_cells2.append(j)

counts=Counter(end_cells2)

print("Count results")
print(counts)

# Select cells present in more than half of the rounds
final_cells=[]
for i in counts.keys():
    if counts[i] >= 6:
        final_cells.append(i)
# If more than one cell of the same type => select only one
df = adata.obs['prediction']
df = df.loc[final_cells]
df = df.drop_duplicates()
final_cells = list(df.index)

pr_res = palantir.core.run_palantir(ms_data, start_cells[0], terminal_states=final_cells, n_jobs = 1)

# Rename trajectories
i=pr_res.branch_probs.columns
df = adata.obs['prediction']
df = df.loc[i]
names = list(df)
pr_res.branch_probs.columns = names

# Plot results
palantir.plot.highlight_cells_on_tsne(umap_seurat, final_cells)
plt.savefig(plot_name + "_umap_final_cells_unsupervised.pdf")

palantir.plot.plot_palantir_results(pr_res, umap_seurat)
plt.savefig(plot_name + "_umap_branches_unsupervised.pdf")

####################################
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
# Impute data
imp_df = palantir.utils.run_magic_imputation(norm_df, dm_res)
imp_df.to_hdf(res_path + "imp_df.hdf5", mode="w", key="imp_df")
