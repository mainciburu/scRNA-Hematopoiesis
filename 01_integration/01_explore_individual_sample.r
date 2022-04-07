###############################################
####### Explore individual samples ############
###############################################

# Input -> 10x matrix
# Create seurat object
# Filter, preprocess, cluster, UMAP
# Plot
###############################################

library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(cluster)
library(future)
source("/home/mainciburu/scRNA/pipeline_2.0/functions.r")

##### INPUT ######
name<-"senior2"
path<-"/home/mainciburu/data/SC_MDS/Senior2/filtered_gene_bc_matrices/GRCh38/"
cols<-"senior2"
condition<-"Elderly"
seurat.file<-"/home/mainciburu/scRNA/senior/seurat_senior2.rds"
########################

# Create seurat object
seurat<-init_seurat(name, path, cols, condition)


# Filter 
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "percent.mito", feature2 = "nFeature_RNA")
VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))
# percent.mito threshold
y<-sapply(seq(1,0,length.out = 41), function(x){sum(seurat$percent.mito<x)})
x<-seq(1,0,length.out = 41)
df<-as.data.frame(cbind(x,y))
ggplot(df, aes(x, y)) + geom_point()

mito.low<-0.005
mito.high<-0.05
count.low<-0
count.high<-35000
feature.low<-500
feature.high<-Inf
seurat<-subset(seurat, percent.mito>mito.low & percent.mito<mito.high & 
			   nCount_RNA > count.low & nCount_RNA < count.high &
			   nFeature_RNA > feature.low & nFeature_RNA < feature.high)

# Activate parallelization
plan("multiprocess", workers = 2)

# Normalize, scale, regress effects and PCA
npcs<-30
vars<-c("S.Score", "G2M.Score")
seurat<-preprocess_seurat(seurat, npcs, vars)

# Decide dimmensionality
ElbowPlot(seurat, ndims = 30)
DimPlot(seurat, reduction = "pca", group.by = "Phase")

# umap
npcs<-15
seurat<-RunUMAP(object = seurat, assay = "RNA", dims = 1:npcs, seed.use = 123)

# explore unwanted effects
DimPlot(seurat, reduction = "umap", group.by = "Phase")
FeaturePlot(seurat, 
            features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
            reduction = "umap", cols = rev(brewer.pal(11, "RdYlBu")))

# clusters
resolution<-0.6
seurat<-cluster_seurat(seurat, npcs, resolution)


# VlnPlots
title<-"Elderly 2"
vlnplot(seurat.obj = seurat, title = title, clusters = "seurat_clusters")

# Cell cycle, counts, features and percent.mito per cluster
VlnPlot(seurat, features = c("S.Score", "G2M.Score", "nCount_RNA", 
                             "nFeature_RNA", "percent.mito"), 
        group.by = "seurat_clusters")

# Lineage markers plots ----------------------------------------------------------

col<-rev(brewer.pal(11, "RdYlBu"))

# B cells
FeaturePlot(seurat, reduction = "umap", 
            features = c("VPREB1", "VPREB3", "CD9", "CD24", "CD79A", "CD79B"), 
            cols = col)

# Dendritic cells
FeaturePlot(seurat, reduction = "umap", 
            features = c("IRF8", "IL3RA", "CLEC4C", "THBD"), 
            cols = col)
# Megakariocytes
FeaturePlot(seurat, reduction = "umap", 
            features = c("PBX1", "HPGDS", "PF4"), 
            cols = col)
# Erythroid
FeaturePlot(seurat, reduction = "umap", 
            features = c("HBB", "HBD", "CA1", "AHSP", "GATA1", "KLF1"), 
            cols = col, pt.size = 0.5)

# HSC
FeaturePlot(seurat, reduction = "umap", 
            features = c("CD34", "CD52", "HOPX", "CD164", "MLLT3", "FOS"), 
            cols = col)

# Granulocytes
FeaturePlot(seurat, reduction = "umap", 
            features = c("CSF3R", "MPO", "SRGN", "CST7"), 
            cols = col)

# Neutrophils
FeaturePlot(seurat, reduction = "umap", 
            features = c("ELANE", "AZU1", "CTSG", "CEBPE"), 
            cols = col)


# Eosinophils
FeaturePlot(seurat, reduction = "umap", 
            features = c("RNASE2", "RNASE3", "PRTN3", "EPX"), 
            cols = col)

# Basophils
FeaturePlot(seurat, reduction = "umap", 
            features = c("HDC", "MS4A2", "MS4A3", "TPSAB1"), 
            cols = col)

# Monocytes
FeaturePlot(seurat, reduction = "umap", 
            features = c("LYZ", "CSTA"), 
            cols = rev(brewer.pal(11, "RdYlBu")))

# Annotate clusters by marker expression
# Labels in scRNA_parameters.xlsx -> cluster_annotation
# young1
#new<-c("LMPP", "GMP_Monocytes", "HSC", "Erythroid", "GMP", "CLP", 
#       "ProB", "MEP", "pDC", "Basophils", "Megakaryocytes")
#new_order<-c("HSC", "LMPP", "GMP", "GMP_Monocytes",
#             "pDC", "CLP", "ProB", "MEP", "Megakaryocytes", "Erythroid", "Basophils")
# young2
#new<-c("LMPP", "GMP1", "GMP_Granulocytes", "Erythroid_early", "CLP_ProB", 
#       "HSC", "pDC", "MEP", "GMP2", "Monocytes", "Erythroid_late", "Basophils")
#new_order<-c("HSC", "LMPP", "GMP_Granulocytes", "GMP1", "GMP2", "Monocytes", "pDC", 
#             "CLP_ProB", "MEP", "Erythroid_early", "Erythroid_late", "Basophils")
# young3
#new<-c("LMPP", "HSC", "GMP_Granulocytes", "pDC", "GMP", "Erythroid", "ProB", "MEP", 
#       "CLP", "Monocytes", "Cycling_LMPP", "Basophils")
#new_order<-c("HSC", "LMPP", "Cycling_LMPP", "GMP_Granulocytes", "GMP", "Monocytes", 
#             "pDC", "CLP", "ProB", "MEP", "Erythroid", "Basophils")

# young4
#new<-c("LMPP", "GMP_Granulocytes", "Erythroid", "HSC", "MEP",
#       "pDC", "GMP", "Cycling_LMPP", "Basophils")
#new_order<-c("HSC", "LMPP", "Cycling_LMPP", "GMP_Granulocytes", "GMP",
#             "pDC", "MEP", "Erythroid", "Basophils")

# young5
#new<-c("LMPP", "Erythroid_early", "GMP_Granulocytes", "MEP", "GMP", "HSC", "ProB",
#       "CLP", "pDC", "Erythroid_late", "T_NK", "Basophils", "Cycling_LMPP", "Megakaryocytes", "Monocytes")
#new_order<-c("HSC", "LMPP", "Cycling_LMPP", "GMP_Granulocytes", "GMP", "Monocytes", 
#             "pDC", "CLP", "ProB", "T_NK", "MEP", "Megakaryocytes", "Erythroid_early", 
#              "Erythroid_late", "Basophils")

# senior1
#new<-c("HSC1", "LMPP", "Erythroid_early", "MEP2", "MEP1",
#       "pDC", "GMP", "Basophils", "Erythroid_late", "HSC2")
#new_order<-c("HSC1", "HSC2", "LMPP", "GMP", "pDC",
#             "MEP1", "MEP2", "Erythroid_early", "Erythroid_late", "Basophils")

# senior2
#new<-c("HSC1", "HSC2", "LMPP", "HSC3", "Erythroid", "MEP", "GMP",
#       "Cycling_LMPP", "pDC", "Basophils")
#new_order<-c("HSC1", "HSC2", "HSC3", "LMPP", "Cycling_LMPP", "GMP",   
#             "pDC", "MEP", "Erythroid", "Basophils")

# senior3
#new<-c("HSC", "LMPP", "Erythroid_early", "MEP", "GMP_Granulocytes", "pDC_Monocytes",
#       "Erythroid_late", "CLP", "Basophils", "Cycling_LMPP1", "ProB", "Cycling_LMPP2",
#       "doublets1?", "doublets2?")
#new_order<-c("HSC", "LMPP", "Cycling_LMPP1", "Cycling_LMPP2", "GMP_Granulocytes",   
#             "pDC_Monocytes", "CLP", "ProB", "MEP", "Erythroid_early", "Erythroid_late",
#             "Basophils", "doublets1?", "doublets2?")


seurat$CellType<-plyr::mapvalues(x = seurat$seurat_clusters, from = levels(seurat$seurat_clusters), to = new)
seurat$CellType<-factor(seurat$CellType, levels = new_order)

vlnplot(seurat, title, "CellType")

# Save object
saveRDS(seurat, file = seurat.file)




