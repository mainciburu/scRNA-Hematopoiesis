##### Integrate samples from the same condition (Young or senior)


library(Seurat)
library(umap)
library(dplyr)
library(RColorBrewer)
library(cluster)

########### FUNCTIONS #######################

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

##############################################

# Load preprocessed seurat objects from individual samples: filtered, normalized, scaled, clustered

seurat1<-readRDS("")
seurat2<-readRDS("")

# List of objects
AllData<-list(seurat1, seurat2)
names(AllData)<-c("seurat1", "seurat2")

# Activate parallelization
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 2097152000)     ## change

# integration
nfeatures<-3000
ndim<-50
int.features<-SelectIntegrationFeatures(object.list = AllData, nfeatures = nfeatures)

anchors<-FindIntegrationAnchors(object.list = AllData, 
                                dims = 1:ndim)
seurat.int<-IntegrateData(anchorset = anchors, dims = 1:ndim)


# Analysis on integrated data

# Recalculate cell Cycle score 
DefaultAssay(seurat.int)<-"RNA"
seurat.int<-CellCycleScoring(seurat.int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
# scale integrated data and regress cc effect
DefaultAssay(seurat.int)<-"integrated"
seurat.int<-ScaleData(seurat.int, vars.to.regress = c("S.Score", "G2M.Score"))

# PCA
seurat.int<-RunPCA(seurat.int, npcs = ndim)
DimPlot(seurat.int, dims = 1:2, reduction = "pca", group.by ="Patient", pt.size = 0.1)
# decide dimensionality
ElbowPlot(seurat.int, ndims = ndim)

# umap
npcs<-25
reduction.name<-"umap"
seurat.int<--RunUMAP(object = seurat.int, assay = "integrated", dims = 1:npcs, reduction.name = reduction.name)
DimPlot(seurat.int, reduction = "umap", dims = 1:2, group.by = "Patient", pt.size = 0.1)

# clustering  !!ONLY WHEN INTEGRATING YOUNG DATA
resolution<-0.4
seurat.int<-cluster_seurat(seurat.int, npcs, resolution)

# VlnPlots
pdf.path<-"/home/mainciburu/scRNA/young/pics/seurat_vln_umap_cluster.pdf"
title<-"Young Integrated"
vlnplot(seurat.int, title, pdf.path)

# save object
saveRDS(seurat.int, file = "seurat_int.rds")
