##### Figure 4 ###########
# UMAP predictions MDS patients
# Proportions barplot
# STREAM plots
# Gene trend heatmaps
# Gene trend lineplots
# GRN heatmap

##########################
library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(plyr)
source("/home/mainciburu/scRNA/colors.r")

### UMAP predictions
mds<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds_integrated.rds")

pdf(file = "/home/mainciburu/scRNA/figures/figure4/UMAP_MDS_celltype.pdf", useDingbats = F,
    width = 12, height = 8)
col<-c(col.young.v3[names(col.young.v3)%in%levels(mds$prediction)], "not assigned"="grey")
DimPlot(mds, reduction = "umap", group.by = "prediction", split.by = "Patient",
               pt.size = 0.5, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()


### Proportions barplot
tt<-prop.table(table(mds$prediction, mds$Patient), margin = 2)
df<-melt(tt)
colnames(df)<-c("Prediction", "Patient", "Proportion")
df$Prediction<-factor(df$Prediction, levels = rev(names(col.young.v3)))
df$Patient<-factor(df$Patient, levels = c("mds10", "mds5", "mds3", "mds1"))
pdf("/home/mainciburu/scRNA/figures/figure4/proportion_mds.pdf", useDingbats = F,
    width = 14, height = 3)
ggplot(df, aes(Proportion, Patient, fill = Prediction)) + geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = rev(col)) + theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold"))
dev.off()


### GSEA senior vs MDS: go to 03_differential_expression/03_GSEA_young_elderly_mds.r

### Gene trends

### SCENIC heatmaps: go to 04_GRN/04_downstream_analysis.r

