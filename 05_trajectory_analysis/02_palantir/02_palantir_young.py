import numpy as np
import pandas as pd     
import pickle
import matplotlib
import matplotlib.pyplot as plt
import scipy
import statsmodels
import scanpy as sc
import palantir
import seaborn as sns
from collections import Counter
import feather


########## Input data ####################
# loom/h5ad object, to extract meta data
h5ad_file="/home/mainciburu/scRNA/palantir/processed/young_int.h5ad"
loom_file="/home/mainciburu/scRNA/palantir/data/young_int.loom"
# Integrated expression matrix (all common genes)
int_df_file="/home/mainciburu/scRNA/palantir/data/young_int_mat.csv"
# Normalized expression matrix (all common genes)
norm_df_full_file="/home/mainciburu/scRNA/palantir/data/young_norm_mat_full.csv"
# UMAP coordinates
umap_seurat_file="/home/mainciburu/scRNA/palantir/data/young_umap.txt"

plot_name = "/home/mainciburu/scRNA/palantir/pics/young"
res_path="/home/mainciburu/scRNA/palantir/results/young/"

##########################################

# load loom // h5ad object 
adata = sc.read_loom(loom_file)
adata.write(filename=h5ad_file)
adata = sc.read(filename=h5ad_file)

# Integrated matrix from csv
int_df=pd.read_csv(int_df_file, index_col=0)
# transpose
int_df=int_df.T

# Import seurat umap coordinates
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

### Select 10 start cells 
# 10 HSC upper on UMAP 2 coordinate
start_cells = ["CAGAGAGCAAGTCTGT_young1", "TGAGGGAGTCACGACC_young5",
"CTCGGAGTCGTCCGTT_young4", "CAAGAAACAAATTGCC_young4",
"TGTTCCGCACATTAGC_young4", "GAGCAGAGTGCAACTT_young4",
"TCGTAGAGTTCCATGA_young4", "GCTACAACATAACGGG_young5",
"TGCTCCAAGCTCCATA_young5", "GTATCTTTCGACCAGC_young4"]

palantir.plot.highlight_cells_on_tsne(umap_seurat, start_cells)
plt.savefig(plot_name + "_umap_start_cells.pdf")

print("Start cells---------------")
for cell in start_cells:
    print(cell)

# Calculate trajectories for each starting cell and store ending cells
ms_data.index = adata.obs.index
end_cells = []
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
df = adata.obs['CellType2']
df = df.loc[final_cells]
df = df.drop_duplicates()
final_cells = list(df.index)

print("final cells")
print(final_cells)

# Repeat trajectories with all ending cells
pr_res = palantir.core.run_palantir(ms_data, start_cells[0], terminal_states=final_cells, n_jobs = 1)

# Rename trajectories
i=pr_res.branch_probs.columns
df = adata.obs['CellType2']
df = df.loc[i]
names = list(df)
pr_res.branch_probs.columns = names

# Plot results
palantir.plot.highlight_cells_on_tsne(umap_seurat, final_cells)
plt.savefig(plot_name + "_umap_final_cells.pdf")

palantir.plot.plot_palantir_results(pr_res, umap_seurat)
plt.savefig(plot_name + "_umap_branches.pdf")

# pseudotime, branch probability, differentiation potential
pr_res.branch_probs.to_csv(res_path + "branch_probs.csv")
pr_res.pseudotime.to_csv(res_path + "pseudotime.csv", )
pr_res.entropy.to_csv(res_path + "diff_potential.csv")

# Save pr_res object
file_pr = open(res_path + 'pr_res.obj', 'wb') 
pickle.dump(pr_res, file_pr)

# Gene expression trends
# !!Use original normalized expression values (not integrated)
norm_df_full = pd.read_csv(norm_df_full_file, index_col = 0)
norm_df_full = pd.DataFrame.transpose(norm_df_full)
# Impute data
imp_df = palantir.utils.run_magic_imputation(norm_df_full, dm_res)
imp_df.to_hdf(res_path + "imp_df.hdf5", mode="w", key="imp_df")
### Gene trends error => continue in R