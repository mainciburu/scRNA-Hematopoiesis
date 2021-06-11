####### Explore individual samples ############
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

############ FUNCTIONS ####################################

# init_seurat - create seurat object and calculate % of mitochondrial genes
# input
	# name: sample name
	# path: path to the filtered 10X matrix
	# cols: string to add to the barcode names
	# condition: Young, Senior or MDS
# output: seurat object 

init_seurat<-function(name, path, cols, condition){
	init.data<-Read10X(data.dir = path)
	colnames(init.data)<-paste0(colnames(init.data), "_", cols)
	seurat.obj<-CreateSeuratObject(counts = init.data, min.cells = 3, min.features = 100, project = name)
	seurat.obj$Patient<-name
	seurat.obj$Condition<-condition

	mito.genes<-grep("MT-", rownames(seurat.obj), value = T)
	percent.mito<-Matrix::colSums(GetAssayData(seurat.obj, "counts")[mito.genes,])/Matrix::colSums(seurat.obj)
	seurat.obj$percent.mito<-percent.mito
	return(seurat.obj)
}

# preprocess_seurat - normalize, calculate cell cycle scores, scale, regress effects, find variable genes and PCA
# input
	# seurat.obj: object
	# npcs: number of principal component to calculate
	# vars: variables from metadata to regress
# output: processed seurat object

preprocess_seurat<-function(seurat.obj, npcs, vars){
	seurat.obj<-NormalizeData(seurat.obj)
	seurat.obj<-CellCycleScoring(object = seurat.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
	seurat.obj<-ScaleData(object = seurat.obj, vars.to.regress = vars)
	seurat.obj<-FindVariableFeatures(seurat.obj, selection.method = "vst",
                            nfeatures = 2000)
	seurat.obj<-RunPCA(seurat.obj, npcs = npcs)
	return(seurat.obj)
}

# cluster_seurat: cluster cells with specific resolution, plot and calculate average silhouette
# input
	# seurat.obj: object
	# npcs: number of significant components
	# resolution: clustering resolution
# output: seurat object with clusters, cluster umap plot and silhouette barplot

cluster_seurat<-function(seurat.obj, npcs, resolution){
	seurat.obj<-FindNeighbors(object = seurat.obj, reduction = "pca", dims = 1:npcs)
	seurat.obj<-FindClusters(seurat.obj, resolution = resolution)
	col<-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.obj$seurat_clusters)))
	print(DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", cols = col))
	# Cluster silhouette
	d<-dist(x = seurat.obj@reductions$pca@cell.embeddings[,1:npcs])
	s<-silhouette(x = as.numeric(seurat.obj$seurat_clusters), dist = d)
	summary(s)
	s.avg<-as.numeric(summary(s)$clus.avg.widths)
	c<-length(unique(seurat.obj$seurat_clusters)) - 1
	barplot(s.avg, horiz = T, names.arg = as.character(0:c), col = col)
	return(seurat.obj)
}

# vlnplot- create umap + vlnplot with lineage markers per cluster
# input
	# seurat.obj: object
	# title: plot title
	# clusters: meta.data column to make groups
# output: plot in pdf

vlnplot<-function(seurat.obj, cluster, title=NULL, col=NULL){
  InfoData<-data.frame(x=seurat.obj@reductions$umap@cell.embeddings[,1], 
                       y=seurat.obj@reductions$umap@cell.embeddings[,2],
                       Cluster=seurat.obj@meta.data[,cluster])
  # Markers
  hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
  lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
  cc<-c("CDC20", "TOP2A")
  gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
  granul<-c("ELANE", "AZU1", "CEBPA", "CEBPE", "CST7")
  mono<-c("LYZ", "CSTA")
  dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4")
  t<-c("JCHAIN", "IKZF1", "CYTH1")
  clp<-c("IL7R", "DNTT")
  prob<-c("VPREB1", "EBF1", "CD79A", "CD79B")
  mep<-c("NFE2", "HFS1", "TAL1")
  mk<-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
  ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
  baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")
  
  markers<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, prob, t, mep, mk, ery, baso)
  markers<-markers[markers%in%rownames(seurat.obj)]
  
  InfoData<-cbind(InfoData,t(as.matrix(seurat.obj@assays$RNA@data[markers,])))
  InfoData2<-melt(data = InfoData, measure.vars = markers)
  if(is.null(col)){
    col<-colorRampPalette(brewer.pal(12, "Paired"))(nrow(unique(seurat.obj[[cluster]])))
  }
  Pvln<-ggplot(InfoData2, aes(x=Cluster, y=value, fill = Cluster))+facet_grid(variable~.)+geom_violin(scale = "width")+ theme(legend.position="none") + labs(x=NULL, y = "Counts")+scale_fill_manual(values=col) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(strip.text.y = element_text(angle = 0)) + theme(strip.background = element_blank()) + scale_y_continuous(limits = c(0,6))
  Pumap<-DimPlot(seurat.obj, reduction = "umap", group.by = cluster, pt.size = 0.5, cols = col) + labs(colour = "Cluster") + ggtitle(label = title)
  plot_grid(Pumap, Pvln, ncol = 1, rel_heights = c(1, 1.5))
}

####################################################


##### INPUT ######
name<-"senior2"
#path<-"/home/mainciburu/data/SC_MDS/Run383/MO258_Count/outs/filtered_feature_bc_matrix/"
path<-"/home/mainciburu/data/SC_MDS/Senior2/filtered_gene_bc_matrices/GRCh38/"
cols<-"senior2"
condition<-"Elderly"
########################

# Create seurat object
seurat<-init_seurat(name, path, cols, condition)


# Filter 
# filtering parameters on scRNA_parameters.xlsx
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
seurat.file<-"/home/mainciburu/scRNA/senior/seurat_senior2.rds"
saveRDS(seurat, file = seurat.file)




