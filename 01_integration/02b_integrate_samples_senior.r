###################################################
######### Integrate elderly samples  ##############
###################################################

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
source("/home/mainciburu/scRNA/pipeline_2.0/functions.r")

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
