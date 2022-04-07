############################
# figure 1
# UMAP unsupervised annotated clusters young 
# UMAP predictions elderly
# Proportions barplot
# Markers DotPlot
# GSEA
############################

library(Seurat)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(reshape)
source("/home/mainciburu/scRNA/colors.r")

young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

# UMAP young CellType2 -------------------------------------------------------------------
pdf(file = "/home/mainciburu/scRNA/figures_dec20/figure1/UMAP_young.pdf", useDingbats = F,
    width = 10, height = 10)
col<-col.young.v3
DimPlot(young, reduction = "umap.int", group.by = "CellType2",
               pt.size = 0.5, cols = col, label = F) + 
  labs(colour = "Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

# UMAP senior predictions ------------------------------------------------------------------
pdf(file = "/home/mainciburu/scRNA/figures_dec20/figure1/UMAP_elderly.pdf", useDingbats = F,
    width = 10, height = 10)
col<-c(col.young.v3[names(col.young.v3)%in%levels(senior$prediction)], "not assigned"="grey")
DimPlot(senior, reduction = "umap.int", group.by = "prediction",
        pt.size = 0.5, cols = col, label = F) + 
  labs(colour = "Predicted Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

# Patient - cluster proportion -------------------------------------------------
# young -> original seurat clusters
tt<-prop.table(table(young$CellType2, young$Patient), margin = 2)
df<-melt(tt)
colnames(df)<-c("CellType2", "Patient", "Proportion")
df$CellType2<-factor(df$CellType2, levels = rev(levels(young$CellType2)))
df$Patient<-factor(df$Patient, levels = c("young5", "young4", "young3", "young2", "young1"))
pdf("/home/mainciburu/scRNA/figures_dec20/figure1/proportion_young.pdf", useDingbats = F,
    width = 14, height = 5)
p<-ggplot(df, aes(Proportion, Patient, fill = CellType2)) + geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = rev(col.young.v3)) + theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold"))
print(p)
dev.off()

# senior -> glmnet predictions
tt<-prop.table(table(senior$prediction, senior$Patient), margin = 2)
df<-melt(tt)
colnames(df)<-c("Prediction", "Patient", "Proportion")
df$Prediction<-factor(df$Prediction, levels = rev(levels(senior$prediction)))
df$Patient<-plyr::mapvalues(x = df$Patient, from = c("senior1", "senior2", "senior3"), 
                            to = c("elderly1", "elderly2", "elderly3"))
df$Patient<-factor(df$Patient, levels = c("elderly3", "elderly2", "elderly1"))
col<-c(col.young.v3[names(col.young.v3)%in%levels(senior$prediction)], "not assigned"="grey")
pdf("/home/mainciburu/scRNA/figures_dec20/figure1/proportion_elderly.pdf", useDingbats = F,
    width = 14, height = 3)
ggplot(df, aes(Proportion, Patient, fill = Prediction)) + geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = rev(col)) + theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold"))
dev.off()

# DotPlot ----------------------------------------------
hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
cc<-c("CDC20", "TOP2A")
gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
granul<-c("ELANE", "AZU1", "CEBPA", "CST7")
mono<-c("LYZ", "CSTA")
dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4", "IGKC", "SPIB", "BLINK")
t<-c("JCHAIN", "PRSS2", "TXNP", "LTB")
nk<-c("CXXC5", "HOXA9", "HOXA10")
clp<-c("IL7R", "DNTT")
prob<-c("VPREB1", "EBF1", "CD79A", "CD79B", "IGHM")
mep<-c("NFE2", "HFS1", "TAL1")
mk<-c("PBX1", "VWF", "FLI1", "ITGA22B", "GP1BA")
ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")

GX<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, prob, t, nk, mep, mk, ery, baso)

GX<-data.frame(Gene=GX,                                                   # Normalized counts for marker genes
               Senior=match(GX,rownames(senior@assays$RNA@data)),
               Young=match(GX,rownames(young@assays$RNA@data)), stringsAsFactors = F)
GX<-GX[-which(is.na(GX$Senior) & is.na(GX$Young)),]

DefaultAssay(young)<-"RNA"
DefaultAssay(senior)<-"RNA"

P1<-DotPlot(young,assay = "RNA",features = rev(GX$Gene),dot.min = 0.4,group.by = "CellType2",
            col.min = -2,col.max = 2,col=c("#f7e69c","red"))

P2<-DotPlot(senior,assay = "RNA",features = rev(GX$Gene),dot.min = 0.4,group.by = "prediction",
            col.min = -2,col.max = 2,col=c("#f7e69c","red"))

YoungData<-P1$data
YoungData$Condition<-"Young"

SeniorData<-P2$data
SeniorData$Condition<-"Senior"

PlotData<-rbind(YoungData,SeniorData)
PlotData$Condition<-factor(PlotData$Condition,levels = c("Young","Senior"))
PlotData$features.plot<-factor(PlotData$features.plot,levels = rev(levels(PlotData$features.plot)))
ggplot(PlotData,aes(x=features.plot,y=id,colour=avg.exp.scaled))+
  geom_point(aes(size=pct.exp))+facet_grid(~Condition)+
  scale_colour_gradient(low = "#f7e69c",high = "red")+theme_classic()+coord_flip()

ggsave(file="/home/mainciburu/scRNA/figures_dec20/figure1/DotPlot_YoungSenior.pdf",
       height = 13, width = 9, useDingbats = F)

# GSEA ----------------------------------------------
# Go to 03_differential_expression/01_GSEA_young_elderly.r
