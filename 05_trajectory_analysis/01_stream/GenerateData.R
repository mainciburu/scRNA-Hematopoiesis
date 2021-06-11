library(Seurat)
source("/home/mainciburu/scRNA/colors.r")
Young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")

data_path = "/home/mainciburu/scRNA/stream/data/"

write.table(colnames(Young@assays$integrated@data),file=paste0(data_path, "barcodes.tsv"),sep="\t",col.names = F,row.names = F,quote = F)
write.table(rownames(Young@assays$integrated@data),file=paste0(data_path, "Genes.tsv"),sep="\t",col.names = F,row.names = F,quote = F)
write.table(Young@assays$integrated@data,file=paste0(data_path, "ExpData.csv"),sep=",",col.names = T,row.names = T)

MetaD<-data.frame(label=Young$CellType2)
MetaD$label_color<-col.young.v3[match(MetaD$label,names(col.young.v3))]
write.table(MetaD,file=paste0(data_path, "metadata.tsv"),sep="\t",col.names =T,row.names = T,quote = F)


Senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

write.table(Senior@assays$integrated@data,file=paste0(data_path, "ExpDataSenior.csv"),sep=",",col.names = T,row.names = T)
write.table(colnames(Senior@assays$integrated@data),file=paste0(data_path, "barcodesSenior.tsv"),sep="\t",col.names = F,row.names = F,quote = F)
write.table(rownames(Senior@assays$integrated@data),file=paste0(data_path, "GeneSenior.tsv"),sep="\t",col.names = F,row.names = F,quote = F)

MetaD<-data.frame(label=Senior$prediction)
MetaD$label_color<-col.young.v3[match(MetaD$label,names(col.young.v3))]
MetaD[which(is.na(MetaD$label_color)),2]<-"#c4c3c2"
write.table(MetaD,file=paste0(data_path, "metadataSenior.tsv"),sep="\t",col.names =T,row.names = T,quote = F)

mds<-readRDS("/home/mainciburu/scRNA/MDS/seurat_mds5.rds")

write.table(mds@assays$RNA@data,file=paste0(data_path, "ExpDataMDS5.csv"),sep=",",col.names = T,row.names = T)
write.table(colnames(mds@assays$RNA@data),file=paste0(data_path, "barcodesMDS5.tsv"),sep="\t",col.names = F,row.names = F,quote = F)
write.table(rownames(mds@assays$RNA@data),file=paste0(data_path, "GeneMDS5.tsv"),sep="\t",col.names = F,row.names = F,quote = F)

MetaD<-data.frame(label=mds$prediction)
MetaD$label_color<-col.young.v3[match(MetaD$label,names(col.young.v3))]
MetaD[which(is.na(MetaD$label_color)),2]<-"#c4c3c2"
write.table(MetaD,file=paste0(data_path, "metadataMDS5.tsv"),sep="\t",col.names =T,row.names = T,quote = F)
