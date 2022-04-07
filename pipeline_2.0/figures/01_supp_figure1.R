############################
# Supplementary figure 1
# UMAP unsupervised annotated clusters individual samples
# UMAP unsupervised annotated clusters integrated by condition
############################

library(Seurat)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
source("/home/mainciburu/scRNA/colors.r")

young1<-readRDS("/home/mainciburu/scRNA/young/seurat_young1.rds")
young2<-readRDS("/home/mainciburu/scRNA/young/seurat_young2.rds")
young3<-readRDS("/home/mainciburu/scRNA/young/seurat_young3.rds")
young4<-readRDS("/home/mainciburu/scRNA/young/seurat_young4.rds")
young5<-readRDS("/home/mainciburu/scRNA/young/seurat_young5.rds")
senior1<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior1.rds")
senior2<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior2.rds")
senior3<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior3.rds")

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_young1.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young1
DimPlot(young1, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_young2.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young2
DimPlot(young2, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_young3.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young3
DimPlot(young3, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_young4.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young4
DimPlot(young4, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_young5.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young5
DimPlot(young5, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_senior1.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.senior1
DimPlot(senior1, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_senior2.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.senior2
DimPlot(senior2, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_senior3.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.senior3
DimPlot(senior3, reduction = "umap", group.by = "CellType",
        pt.size = 0.8, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

rm(young1, young2, young3, young4, young5, senior1, senior2, senior3)

##### integrated umaps
young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")
pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_youngIntegrated.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young.v3
DimPlot(young, reduction = "umap.int", group.by = "CellType2",
        pt.size = 0.5, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure1/UMAP_seniorIntegrated.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.senior.v4
DimPlot(senior, reduction = "umap.int", group.by = "CellType2",
        pt.size = 0.5, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()
