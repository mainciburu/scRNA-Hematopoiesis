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
source("/home/mainciburu/scRNA/colors.r")
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
seurat1<-readRDS("/home/mainciburu/scRNA/young/seurat_young1.rds")
seurat2<-readRDS("/home/mainciburu/scRNA/young/seurat_young2.rds")
seurat3<-readRDS("/home/mainciburu/scRNA/young/seurat_young3.rds")
seurat4<-readRDS("/home/mainciburu/scRNA/young/seurat_young4.rds")
seurat5<-readRDS("/home/mainciburu/scRNA/young/seurat_young5.rds")
# List of objects
AllData<-list(young1=seurat1, young2=seurat2, young3=seurat3, young4=seurat4, young5=seurat5)

# Activate parallelization
plan("multiprocess", workers = 8)
options(future.globals.maxSize= 6000000000)     ## change

############ integration ####################
ndim<-50

anchors<-FindIntegrationAnchors(object.list = AllData, 
                                dims = 1:ndim
                                )
####### integrate common genes in young and senior ############
young1<-readRDS("/home/mainciburu/scRNA/young/seurat_young1.rds")
young2<-readRDS("/home/mainciburu/scRNA/young/seurat_young2.rds")
young3<-readRDS("/home/mainciburu/scRNA/young/seurat_young3.rds")
young4<-readRDS("/home/mainciburu/scRNA/young/seurat_young4.rds")
young5<-readRDS("/home/mainciburu/scRNA/young/seurat_young5.rds")
senior1<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior1.rds")
senior2<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior2.rds")
g<-intersect(x = intersect(rownames(young1), rownames(young2)),
             y = intersect(rownames(young3), rownames(young4)))
g<-intersect(x = intersect(rownames(senior1), rownames(senior2)),
             y = g)
write.table(g, file = "/home/mainciburu/scRNA/common_genes.txt", 
			  quote = FALSE, row.names = FALSE, col.names=FALSE)
rm(young1, young2, young3, young4, senior1, senior2)
#############################

g<-read.table("/home/mainciburu/scRNA/common_genes.txt", header = F)
g<-as.character(g$V1)
seurat.int<-IntegrateData(anchorset = anchors, dims = 1:ndim, features.to.integrate = g)
#############################

########## Analysis on integrated data ##############
# Recalculate cell Cycle score 
DefaultAssay(seurat.int)<-"RNA"
seurat.int<-CellCycleScoring(seurat.int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# scale integrated data and regress cc effect
DefaultAssay(seurat.int)<-"integrated"
seurat.int<-ScaleData(seurat.int, vars.to.regress = c("S.Score", "G2M.Score"), features = g)

# PCA
seurat.int<-RunPCA(seurat.int, npcs = ndim)
DimPlot(seurat.int, dims = 1:2, reduction = "pca", group.by ="Patient", pt.size = 0.1)
# decide dimensionality
ElbowPlot(seurat.int, ndims = ndim)

# umap
npcs<-35
reduction.name<-"umap"
seurat.int<-RunUMAP(object = seurat.int, assay = "integrated", dims = 1:npcs, reduction.name = reduction.name, seed.use = 123)
DimPlot(seurat.int, reduction = "umap", dims = 1:2, group.by = "Patient", pt.size = 0.1)

# clustering  
resolution<-0.5
seurat.int<-cluster_seurat(seurat.int, npcs, resolution)

# VlnPlots
title<-"Young Integrated"
pdf("young.pdf", width = 10, height = 15, useDingbats = F)
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

#nCounts and nFeatures per cluster
VlnPlot(seurat.int, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "G2M.Score", "S.Score"), 
        group.by = "seurat_clusters")

######## Transfer from greenleaf data #############
gr<-readRDS("/home/mainciburu/scRNA/greenleaf/seurat_greenleaf_cd34.rds")
gr<-ScaleData(gr, features = rownames(gr))
gr<-FindVariableFeatures(gr, selection.method = "vst")
gr<-RunPCA(gr, npcs = 50)
nfeatures<-3000
ndim<-30
DefaultAssay(seurat.int)<-"RNA"
transfer.features<-SelectIntegrationFeatures(object.list = list(gr, seurat.int), nfeatures = nfeatures)
anchors<-FindTransferAnchors(reference = gr, query = seurat.int, 
                             dims = 1:ndim, reference.assay = "RNA", query.assay = "RNA", 
                             features = transfer.features)
predictions <- TransferData(anchorset = anchors, refdata = gr$BioClassification, 
                            dims = 1:ndim)
seurat.int$greenleaf<-predictions$predicted.id
seurat.int$greenleaf[predictions$prediction.score.max<0.5]<-"not assigned"

df2<-sapply(levels(seurat.int$seurat_clusters), function(i){
  maxlab<-names(which.max(table(seurat.int$greenleaf[seurat.int$seurat_clusters==i])))
  pct<-sum(seurat.int$greenleaf[seurat.int$seurat_clusters==i]==maxlab)/sum(seurat.int$seurat_clusters==i)
  return(c(maxlab, pct))
}
)
df<-rbind(df, df2)
t(df)

p1<-DimPlot(seurat.int, reduction = "umap", group.by = "greenleaf", cols = col.greenleaf)
p2<-DimPlot(seurat.int, reduction = "umap", group.by = "seurat_clusters", 
            cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.int$seurat_clusters))))
cowplot::plot_grid(p1, p2, nrow = 1)
#################################################

# Cluster markers
markers<-FindAllMarkers(seurat.int, assay = "RNA", features = g, only.pos = T, test.use = "MAST")

VlnPlot(seurat.int, features = c("G2M.Score", "S.Score", "nCount_RNA", "nFeature_RNA", "percent.mito"),
        group.by = "CellType2", pt.size = 0, cols = col.young.v3)

# Cell cycle phases
tt<-prop.table(table(seurat.int$CellType2, seurat.int$Phase), margin = 1)
df<-melt(tt)
df$CellType<-factor(df$CellType, levels = levels(seurat.int$CellType2))
colnames(df)<-c("CellType", "Phase", "Proportion")
ggplot(df, aes(CellType, Proportion, fill = Phase)) + geom_bar(stat = "identity", position = position_dodge())


# Name integrated clusters
new.labs<-c("LMPP", "Erythroid_early", "MEP", "GMP", "GMP_Granulocytes",
            "HSC", "pDC1", "ProB", "Cycling_LMPP1", "Erythroid_late", "Basophils",
            "pDC2", "CLP", "Cycling_LMPP2", "Monocytes", "T_NK", "Megakaryocytes",
            "ProB2", "Erythroid_late2", "Monocytes2")
seurat.int$CellType2.orig<-plyr::mapvalues(x = seurat.int$seurat_clusters, 
                                from = levels(seurat.int$seurat_clusters),
                                to = new.labs)
seurat.int$CellType2.orig<-factor(seurat.int$CellType2.orig, 
                             levels = c("HSC", "LMPP", "Cycling_LMPP1", "Cycling_LMPP2", "GMP",
                                        "GMP_Granulocytes", "Monocytes", "Monocytes2", "pDC1", "pDC2",
                                        "CLP", "ProB", "ProB2", "T_NK","MEP", "Megakaryocytes", "Erythroid_early",
                                        "Erythroid_late", "Erythroid_late2", "Basophils"))

# VlnPlot
pdf("scRNA/young/pics/young_v3/umap_vln_orig.pdf", width = 10, height = 25, useDingbats = F)
vlnplot(seurat.obj = seurat.int, cluster = "CellType2.orig", title = "Young", col = col.young.v3.orig)
dev.off()

# Clusters 17:19 (ProB2, Erythroid_late2, Monocytes2)
# Small clusters, very similar to the one next to them in umap, but with lower counts/features
VlnPlot(seurat.int, features = c("G2M.Score", "S.Score", "nCount_RNA", "nFeature_RNA", "percent.mito"),
        group.by = "seurat_clusters", pt.size = 0)

# Remove
Idents(seurat.int)<-"seurat_clusters"
seurat.int<-subset(seurat.int, idents=0:16)

# Merge similar clusters
# LMPPs (LMPP, cycling LMPP1, cycling LMPP2)
seurat.int$CellType2<-seurat.int$CellType2.orig
seurat.int$CellType2<-as.character(seurat.int$CellType2)
seurat.int$CellType2[seurat.int$CellType2=="Cycling_LMPP1"]<-"LMPP"
seurat.int$CellType2[seurat.int$CellType2=="Cycling_LMPP2"]<-"LMPP"

# pDC (pDC1 and pDC2)
seurat.int$CellType2[seurat.int$CellType2=="pDC1"]<-"pDC"
seurat.int$CellType2[seurat.int$CellType2=="pDC2"]<-"pDC"

seurat.int$CellType2<-factor(seurat.int$CellType2, levels = c("HSC", "LMPP", "GMP", "GMP_Granulocytes", "Monocytes", "pDC",
                                                              "CLP", "ProB", "T_NK","MEP", "Megakaryocytes", "Erythroid_early",
                                                              "Erythroid_late", "Basophils"))

# VlnPlot
pdf("scRNA/young/pics/young_v3/umap_vln_merged.pdf", width = 10, height = 25, useDingbats = F)
vlnplot(seurat.obj = seurat.int, cluster = "CellType2", title = "Young", col = col.young.v3)
dev.off()

# Remove meta.data columns with NA
seurat.int@meta.data<-seurat.int@meta.data[,colSums(is.na(seurat.int@meta.data))==0]

# Save
seurat.file<-"/home/mainciburu/scRNA/young/seurat_young_v3.rds"
saveRDS(object = seurat.int, file = seurat.file)
