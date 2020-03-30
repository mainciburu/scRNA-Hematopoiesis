
library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(umap)
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

# preprocess_seurat - normalize, calculate cell cycle scores, scale, regress effects and PCA
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
	DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", cols = col)
	# Cluster silhouette
	d<-dist(x = seurat.obj@reductions$pca@cell.embeddings[,1:npcs])
	s<-silhouette(x = as.numeric(seurat.obj$seurat_clusters), dist = d)
	summary(s)
	s.avg<-as.numeric(summary(s)$clus.avg.widths)
	c<-length(unique(seurat.obj$seurat_clusters))
	barplot(s.avg, horiz = T, names.arg = as.character(0:c), col = col)
	return(seurat.obj)
}

# vlnplot- create umap + vlnplot with lineage markers per cluster
# input
	# seurat.obj: object
	# title: plot title
	# pdf.path: path and name for the pdf
# output: plot in pdf

vlnplot<-function(seurat.obj, title, pdf.path){
	InfoData<-data.frame(x=seurat.obj@reductions$umap@cell.embeddings[,1], 
                     y=seurat.obj@reductions$umap@cell.embeddings[,2],
                     Cluster=seurat.obj@meta.data$seurat_clusters)
	# Markers
	hsc<-c("CRHBP", "HOPX")
	lmpp<-c("CD34")  
	cc<-c("TOP2A", "CDC20")
	gmp<-c("MPO", "CTSG", "PRTN3")
	granul<-c("ELANE")
	mono<-c("LYZ", "CSTA")
	dc<-c("IL3RA", "IRF8")
	clp<-c("IL7R", "DNTT")
	preb<-c("VPREB1", "EBF1")
	mep<-c("PBX1", "VWF")
	ery<-c("HBD", "CA1")
	baso<-c("HDC", "MS4A2")
	markers<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, preb, cmp, mep, ery, baso)
	markers<-markers[markers%in%rownames(seurat.obj)]

	InfoData<-cbind(InfoData,t(as.matrix(seurat.obj@assays$RNA@data[markers,])))
	InfoData2<-melt(data = InfoData, measure.vars = markers)

	col<-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.obj$seurat_clusters)))
	Pvln<-ggplot(InfoData2, aes(x=Cluster, y=value, fill = Cluster))+facet_grid(variable~.)+geom_violin(scale = "width")+ theme(legend.position="none") + labs(x=NULL, y = "Counts")+scale_fill_manual(values=col) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(strip.text.y = element_text(angle = 0)) + theme(strip.background = element_blank()) + scale_y_continuous(limits = c(0,6))
	Pumap<-DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, cols = col) + labs(colour = "Cluster") + ggtitle(label = title)
	pdf(pdf.path, height = 15, useDingbats = F)
	plot_grid(Pumap, Pvln, ncol = 1)
	dev.off()

}


####################################################


##### INPUT ######
name<-"mo258"
path<-"data/SC_MDS/Run383/mo258_Count/outs/filtered_feature_bc_matrix/"
cols<-"mo258"
condition<-"Young"
########################

# Create seurat object
seurat<-init_seurat(name, path, cols, condition)


# Filter 
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "percent.mito", feature2 = "nFeature_RNA")
VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))

mito.low<-0.005
mito.high<-0.05
count.low<-0
count.high<-30000
feature.low<-0
feature.high<-Inf
seurat<-subset(seurat, percent.mito>mito.low & percent.mito<mito.high & 
			   nCount_RNA > count.low & nCount_RNA < count.high &
			   nFeature_RNA > feature.low & nFeature_RNA < feature.high)

# Activate parallelization
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 2097152000)

# Normalize, scale, regress effects and PCA
npcs<-30
vars<-c("S.Score", "G2M.Score")
seurat<-preprocess_seurat(seurat, npcs, vars)

# Decide dimmensionality
ElbowPlot(seurat, ndims = 30)
DimPlot(seurat, reduction = "pca", group.by = "Phase")

# umap
npcs<-15
seurat<-RunUMAP(object = seurat, assay = "RNA", dims = 1:npcs)

# explore unwanted effects
DimPlot(seurat, reduction = "umap", group.by = "Phase")
FeaturePlot(seurat, 
            features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
            reduction = "umap", cols = rev(brewer.pal(11, "RdYlBu")))

# clusters
resolution<-0.6
seurat<-cluster_seurat(seurat, npcs, resolution)


# VlnPlots
pdf.path<-"/home/mainciburu/scRNA/young/pics/seurat_vln_umap_cluster.pdf"
title<-"Young 4"
vlnplot(seurat, title, pdf.path)


# Save object
saveRDS(seurat, file = "seurat_seurat.rds")