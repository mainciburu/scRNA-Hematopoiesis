######### Figure 3 #################
# SCENIC heatmaps
# UMAP expression and AUC
# SCENIC networks
# Enrichment barplots
#########################################

library(Seurat)
library(RColorBrewer)
library(ggsignif)
library(ggplot2)
source("/home/mainciburu/scRNA/colors.r")

################
young<-readRDS("scRNA/young/seurat_young_v3.rds")
senior<-readRDS("scRNA/senior/seurat_senior_v4.rds")





