##############################################
######### Integrate MDS samples ##############
##############################################

# Input -> preprocessed seurat objects
# Integrate, rescale, PCA, cluster, UMAP
# Plot
##############################################

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
source("/home/mainciburu/scRNA/pipeline_2.0/functions.r")

# Load preprocessed seurat objects from individual samples: filtered, normalized, scaled, clustered
seurat1<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds1.rds")
seurat2<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds2.rds")
seurat3<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds3.rds")
seurat4<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds4.rds")

plot_path<-"/home/mainciburu/scRNA/MDS_paper/plots/integration/"
seurat.file<-"/home/mainciburu/scRNA/MDS_paper/seurat_mds_integrated.rds"

AllData<-list(mds1=seurat1, mds2=seurat2, mds3=seurat3, mds4=seurat4)
rm(seurat1, seurat2, seurat3, seurat4)

# Activate parallelization
plan("multiprocess", workers = 8)
options(future.globals.maxSize= 6000000000)     ## change

############ integration ####################
ndim<-50

anchors<-FindIntegrationAnchors(object.list = AllData, 
                                dims = 1:ndim
                                )
# integrate common genes 
g<-read.table("/home/mainciburu/scRNA/common_genes.txt", header = F)
g<-as.character(g$V1)
g<-intersect(g, rownames(AllData[[1]]))
g<-intersect(g, rownames(AllData[[2]]))
g<-intersect(g, rownames(AllData[[3]]))
g<-intersect(g, rownames(AllData[[4]]))
rm(AllData)

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
pdf(paste0(plot_path, "pca_patient.pdf"))
DimPlot(seurat.int, dims = 1:2, reduction = "pca", group.by ="Patient", pt.size = 0.1)
dev.off()
# decide dimensionality
pdf(paste0(plot_path, "elbow_plot.pdf"))
ElbowPlot(seurat.int, ndims = ndim)
dev.off()

# umap
npcs<-30
reduction.name<-"umap"
seurat.int<-RunUMAP(object = seurat.int, assay = "integrated", dims = 1:npcs, reduction.name = reduction.name, seed.use = 123)
pdf(paste0(plot_path, "umap_patient.pdf"))
DimPlot(seurat.int, reduction = "umap", dims = 1:2, group.by = "Patient", pt.size = 0.1)
dev.off()
# clustering  
#resolution<-0.5
#seurat.int<-cluster_seurat(seurat.int, npcs, resolution)

# Remove meta.data columns with NA
seurat.int@meta.data<-seurat.int@meta.data[,colSums(is.na(seurat.int@meta.data))==0]

# Factor predicted celltypes
seurat.int$prediction<-factor(seurat.int$prediction, levels = c("HSC", "LMPP", "GMP", "GMP_Granulocytes", "Monocytes", "pDC",
                                                              "CLP", "ProB", "T_NK","MEP", "Megakaryocytes", "Erythroid_early",
                                                              "Erythroid_late", "Basophils", "not assigned"))

# Save
saveRDS(object = seurat.int, file = seurat.file)

# VlnPlots
title<-"MDS Integrated"
pdf(paste0(plot_path, "prediction_markers_vln_umap.pdf"), width = 10, height = 15, useDingbats = F)
vlnplot(seurat.obj = seurat.int, title = title, cluster = "prediction", col = col.young.v3)
dev.off()

#nCounts and nFeatures per cluster
pdf(paste0(plot_path, "unwanted_effects_vln.pdf"), width = 10, height = 15, useDingbats = F)
VlnPlot(seurat.int, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "G2M.Score", "S.Score"), 
        group.by = "prediction")
dev.off()

# UMAP separated per patient
pdf(paste0(plot_path, "umap_patient_celltype.pdf"))
DimPlot(seurat.int, reduction = "umap", dims = 1:2, group.by = "prediction", split.by="Patient", pt.size = 0.1)
dev.off()


# Save UMAP per patient
umap<-data.frame(seurat.int@reductions$umap@cell.embeddings)
umap<-cbind(umap, Patient=seurat.int$Patient)
rm(seurat.int)

umap.x<-umap[umap$Patient=="mds1",1:2]
seurat<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds1.rds")
seurat[["umap.int"]] <- CreateDimReducObject(embeddings = as.matrix(umap.x), key = "UMAP_", assay = DefaultAssay(seurat))
saveRDS(seurat, file = "/home/mainciburu/scRNA/MDS_paper/seurat_mds1.rds")
write.table(umap.x, file = "/home/mainciburu/scRNA/MDS_paper/mds1_umap_integrated.txt")

umap.x<-umap[umap$Patient=="mds2",1:2]
seurat<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds2.rds")
seurat[["umap.int"]] <- CreateDimReducObject(embeddings = as.matrix(umap.x), key = "UMAP_", assay = DefaultAssay(seurat))
saveRDS(seurat, file = "/home/mainciburu/scRNA/MDS_paper/seurat_mds2.rds")
write.table(umap.x, file = "/home/mainciburu/scRNA/MDS_paper/mds2_umap_integrated.txt")

umap.x<-umap[umap$Patient=="mds3",1:2]
seurat<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds3.rds")
seurat[["umap.int"]] <- CreateDimReducObject(embeddings = as.matrix(umap.x), key = "UMAP_", assay = DefaultAssay(seurat))
saveRDS(seurat, file = "/home/mainciburu/scRNA/MDS_paper/seurat_mds3.rds")
write.table(umap.x, file = "/home/mainciburu/scRNA/MDS_paper/mds3_umap_integrated.txt")

umap.x<-umap[umap$Patient=="mds4",1:2]
seurat<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds4.rds")
seurat[["umap.int"]] <- CreateDimReducObject(embeddings = as.matrix(umap.x), key = "UMAP_", assay = DefaultAssay(seurat))
saveRDS(seurat, file = "/home/mainciburu/scRNA/MDS_paper/seurat_mds4.rds")
write.table(umap.x, file = "/home/mainciburu/scRNA/MDS_paper/mds4_umap_integrated.txt")
