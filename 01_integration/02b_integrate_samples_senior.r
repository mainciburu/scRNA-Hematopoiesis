######### Integrate samples from the same condition ##############
# Input -> preprocessed seurat objects
# Integrate, rescale, PCA, cluster, UMAP
# Plot
#################################################################

library(Seurat)
library(cowplot)
library(dplyr)
library(RColorBrewer)
library(cluster)
library(future)
library(ggplot2)
library(reshape)
library(plyr)
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
  mono<-c("LYZ", "CSTA", "CD14")
  dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4", "ADA", "TCF4", "IGKC", "SPIB", "PLAC8", "LILRA4", "GZMB", "BLINK")
  t<-c("JCHAIN", "IKZF1", "CYTH1", "PRSS2", "CD52", "TXNP", "LTB")
  nk<-c("TSC22D1", "CXXC5", "HOXA9", "HOXA10")
  clp<-c("IL7R", "DNTT")
  prob<-c("VPREB1", "EBF1", "CD79A", "CD79B", "IGHM")
  mep<-c("NFE2", "HFS1", "TAL1")
  mk<-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
  ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
  baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")
  
  markers<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, prob, t, nk, mep, mk, ery, baso)
  markers<-markers[markers%in%rownames(seurat.obj)]
  
  InfoData<-cbind(InfoData,t(as.matrix(seurat.obj@assays$RNA@data[markers,])))
  InfoData2<-melt(data = InfoData, measure.vars = markers)
  if(is.null(col)){
    col<-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.obj$seurat_clusters)))
  }
  Pvln<-ggplot(InfoData2, aes(x=Cluster, y=value, fill = Cluster))+facet_grid(variable~.)+geom_violin(scale = "width")+ theme(legend.position="none") + labs(x=NULL, y = "Counts")+scale_fill_manual(values=col) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(strip.text.y = element_text(angle = 0)) + theme(strip.background = element_blank()) + scale_y_continuous(limits = c(0,6))
  Pumap<-DimPlot(seurat.obj, reduction = "umap", group.by = cluster, pt.size = 0.5, cols = col) + labs(colour = "Cluster") + ggtitle(label = title)
  plot_grid(Pumap, Pvln, ncol = 1, rel_heights = c(1, 1.5))
}


##############################################

# Load preprocessed seurat objects from individual samples: filtered, normalized, scaled, clustered

seurat1<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior1.rds")
seurat2<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior2.rds")
seurat3<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior3.rds")

# List of objects
AllData<-list(senior1=seurat1, senior2=seurat2, senior3=seurat3)
rm(seurat1, seurat2, seurat3)

# Activate parallelization
plan("multiprocess", workers = 8)
options(future.globals.maxSize= 10000000000)     ## change

# integration
ndim<-50

anchors<-FindIntegrationAnchors(object.list = AllData, 
                                dims = 1:ndim
                                )

# integrate every common gene in young and senior
g<-read.table("/home/mainciburu/scRNA/common_genes.txt", header = F)
g<-as.character(g$V1)
seurat.int<-IntegrateData(anchorset = anchors, dims = 1:ndim, features.to.integrate = g)


# Analysis on integrated data

# Recalculate cell Cycle score 
DefaultAssay(seurat.int)<-"RNA"
seurat.int<-CellCycleScoring(seurat.int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# scale integrated data and regress cc effect
DefaultAssay(seurat.int)<-"integrated"
#plan("multiprocess", workers = 6)
seurat.int<-ScaleData(seurat.int, vars.to.regress = c("S.Score", "G2M.Score"), features = g)

# PCA
seurat.int<-RunPCA(seurat.int, npcs = ndim)
DimPlot(seurat.int, dims = 1:2, reduction = "pca", group.by ="Patient", pt.size = 0.1)
# decide dimensionality
ElbowPlot(seurat.int, ndims = ndim)

# umap
npcs<-21
reduction.name<-"umap"
seurat.int<-RunUMAP(object = seurat.int, assay = "integrated", dims = 1:npcs, 
                    reduction.name = reduction.name, seed.use = 123)
DimPlot(seurat.int, reduction = "umap", dims = 1:2, group.by = "Patient", pt.size = 0.1)

# clustering  
resolution<-0.3
seurat.int<-cluster_seurat(seurat.int, npcs, resolution)

# VlnPlots
title<-"Senior Integrated"
pdf("senior_integrated.pdf", width = 8, height = 16)
vlnplot(seurat.obj = seurat.int, title = title, cluster = "seurat_clusters")
dev.off()

# individual vs integrated clusters
df<-sapply(levels(seurat.int$seurat_clusters), function(i){
  maxlab<-names(which.max(table(seurat.int$CellType[seurat.int$seurat_clusters==i])))
  pct<-sum(seurat.int$CellType[seurat.int$seurat_clusters==i]==maxlab)/sum(seurat.int$seurat_clusters==i)
  return(c(maxlab, pct))
}
)
t(df)


# Name clusters
new.labs<-c("HSC", "Erythroid", "LMPP", "MEP",  
            "pDC_Monocytes", "GMP_Granulocytes", "Cycling_LMPP", "Basophils", "ProB")
seurat.int$CellType2<-plyr::mapvalues(x = seurat.int$seurat_clusters, 
                                from = levels(seurat.int$seurat_clusters),
                                to = new.labs)
seurat.int$CellType2<-factor(seurat.int$CellType2, 
                             levels = c("HSC", "LMPP", "Cycling_LMPP", "GMP_Granulocytes", 
                                        "pDC_Monocytes", "ProB", "MEP", "Erythroid",
                                        "Basophils"))

# VlnPlot
title<-"Senior Integrated"
pdf("senior_integrated.pdf", width = 8, height = 16)
vlnplot(seurat.obj = seurat.int, cluster = "CellType2", title = title)
dev.off()

# save 
seurat.file<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"
saveRDS(object = seurat.int, file = seurat.file)
