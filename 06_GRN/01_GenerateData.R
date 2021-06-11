library(Seurat)
source("/home/mainciburu/scRNA/colors.r")

# Young
young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")

g<-read.table("/home/mainciburu/scRNA/common_genes.txt", stringsAsFactors=F)$V1

x<-young@assays$RNA@data[g,]
write.csv(x = x, file = "/home/mainciburu/scRNA/scenic/data/young_norm_mat_full.csv")

MetaD<-data.frame(CellType=young$CellType2)
write.table(MetaD,file="/home/mainciburu/scRNA/scenic/data/young_CellType.txt",col.names =T,row.names = T)


# Senior
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

x<-senior@assays$RNA@data[g,]
write.csv(x = x, file = "/home/mainciburu/scRNA/scenic/data/senior_norm_mat_full.csv")

MetaD<-data.frame(CellType=senior$prediction)
write.table(MetaD,file="/home/mainciburu/scRNA/scenic/data/senior_CellType.txt",col.names =T,row.names = T)

# AML
aml<-readRDS("/home/mainciburu/scRNA/MDS/seurat_aml2.rds")

g<-g[g%in%rownames(aml)]

x<-aml@assays$RNA@data[g,]
write.csv(x = x, file = "/home/mainciburu/scRNA/scenic/data/aml_norm_mat_full.csv")

MetaD<-data.frame(CellType=aml$prediction)
write.table(MetaD,file="/home/mainciburu/scRNA/scenic/data/aml_CellType.txt",col.names =T,row.names = T)
