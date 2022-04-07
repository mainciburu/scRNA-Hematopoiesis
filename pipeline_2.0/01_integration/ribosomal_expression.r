
library(Seurat)
library(umap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

seurat1.file<-"/home/mainciburu/scRNA/young/seurat_young_v3.rds"
seurat2.file<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"

young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

DefaultAssay(young)<-"RNA"
DefaultAssay(senior)<-"RNA"

young[["percent.rp"]] <- PercentageFeatureSet(young, pattern = "^RPS") + PercentageFeatureSet(young, pattern = "^RPL")
senior[["percent.rp"]] <- PercentageFeatureSet(senior, pattern = "^RPS") + PercentageFeatureSet(senior, pattern = "^RPL")

df.y<-data.frame(Condition = young$Condition,
				 Patient = young$Patient,
				 RibosomalExpression = young$percent.rp,
				 UMAP1 = young@reductions$umap.int@cell.embeddings[,1],
				 UMAP2 = young@reductions$umap.int@cell.embeddings[,2],
				 CellType = young$CellType2,
				 row.names = colnames(young))
rm(young)
df.s<-data.frame(Condition = senior$Condition,
				 Patient = senior$Patient,
				 RibosomalExpression = senior$percent.rp,
				 UMAP1 = senior@reductions$umap.int@cell.embeddings[,1],
				 UMAP2 = senior@reductions$umap.int@cell.embeddings[,2],
				 CellType = senior$prediction,
				 row.names = colnames(senior))
rm(senior)

mds1<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds1.rds")
mds2<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds3.rds")
mds3<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds5.rds")
mds4<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds10.rds")

mds1[["percent.rp"]] <- PercentageFeatureSet(mds1, pattern = "^RPS") + PercentageFeatureSet(mds1, pattern = "^RPL")
mds2[["percent.rp"]] <- PercentageFeatureSet(mds2, pattern = "^RPS") + PercentageFeatureSet(mds2, pattern = "^RPL")
mds3[["percent.rp"]] <- PercentageFeatureSet(mds3, pattern = "^RPS") + PercentageFeatureSet(mds3, pattern = "^RPL")
mds4[["percent.rp"]] <- PercentageFeatureSet(mds4, pattern = "^RPS") + PercentageFeatureSet(mds4, pattern = "^RPL")

df.m1<-data.frame(Condition = mds1$Condition,
				 Patient = mds1$Patient,
				 RibosomalExpression = mds1$percent.rp,
				 UMAP1 = mds1@reductions$umap.int@cell.embeddings[,1],
				 UMAP2 = mds1@reductions$umap.int@cell.embeddings[,2],
				 CellType = mds1$prediction,
				 row.names = colnames(mds1))

df.m2<-data.frame(Condition = mds2$Condition,
				 Patient = mds2$Patient,
				 RibosomalExpression = mds2$percent.rp,
				 UMAP1 = mds2@reductions$umap.int@cell.embeddings[,1],
				 UMAP2 = mds2@reductions$umap.int@cell.embeddings[,2],
				 CellType = mds2$prediction,
				 row.names = colnames(mds2))

df.m3<-data.frame(Condition = mds3$Condition,
				 Patient = mds3$Patient,
				 RibosomalExpression = mds3$percent.rp,
				 UMAP1 = mds3@reductions$umap.int@cell.embeddings[,1],
				 UMAP2 = mds3@reductions$umap.int@cell.embeddings[,2],
				 CellType = mds3$prediction,
				 row.names = colnames(mds3))

df.m4<-data.frame(Condition = mds4$Condition,
				 Patient = mds4$Patient,
				 RibosomalExpression = mds4$percent.rp,
				 UMAP1 = mds4@reductions$umap.int@cell.embeddings[,1],
				 UMAP2 = mds4@reductions$umap.int@cell.embeddings[,2],
				 CellType = mds4$prediction,
				 row.names = colnames(mds4))

rm(mds1,mds2,mds3,mds4)

df<-rbind(df.y,df.s,df.m1,df.m2,df.m3,df.m4)

### UMAP

pdf("scRNA/MDS_paper/plots/ribosomal_expression.pdf", height = 12, width = 10)
ggplot(df, aes(UMAP1, UMAP2, color = RibosomalExpression)) + geom_point(size = 0.3) + facet_wrap(~Patient, nrow = 4) + 
		theme_bw() + scale_color_gradient2(low = "grey", high = "red") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


### Violin per celltype
pdf("scRNA/MDS_paper/plots/ribosomal_expression_vln.pdf", height = 10, width = 10)
ggplot(df, aes(Patient, RibosomalExpression, fill = CellType)) + geom_violin(scale = "width") + facet_wrap(~CellType, nrow = 4) + 
		theme_bw() + scale_fill_manual(values = col.young.v3) + theme()
dev.off()

