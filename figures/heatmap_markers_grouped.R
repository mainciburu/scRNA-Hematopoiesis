# Heatmap
hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
cc<-c("CDC20", "TOP2A")
gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
granul<-c("ELANE", "AZU1", "CEBPA", "SRGN", "CEBPE", "CST7")
mono<-c("LYZ", "CSTA")
dc<-c("IRF8", "IRF7", "IL3RA", "CD74", "CLEC4")
t<-c("JCHAIN")
clp<-c("IL7R", "DNTT")
prob<-c("VPREB1", "EBF1", "CD79A", "CD79B")
mep<-c("NFE2", "HFS1", "TAL1")
mk<-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")

markers<-c(hsc, lmpp, cc, gmp, granul, mono, dc, t, clp, prob, mep, mk, ery, baso)

#young<-ScaleData(young, assay = "RNA")
#senior<-ScaleData(senior, assay = "RNA")

markers<-markers[markers%in%rownames(young@assays$RNA@scale.data)]
markers<-markers[markers%in%rownames(senior@assays$RNA@scale.data)]

infoData<-data.frame(cluster=young$CellType2.orig, row.names = rownames(young@meta.data))
infoData<-cbind(infoData, t(young@assays$RNA@scale.data[markers,]))
for(g in markers){
  infoData[[g]]<-ave(infoData[[g]], infoData$cluster, FUN = mean)
}
meanData<-infoData[!duplicated(infoData$CRHBP),]
mat<-t(meanData[,2:ncol(meanData)])
colnames(mat)<-as.character(unique(infoData$cluster))
mat<-mat[,levels(infoData$cluster)]
mat1<-mat

infoData<-data.frame(cluster=senior$CellType2, row.names = rownames(senior@meta.data))
infoData<-cbind(infoData, t(senior@assays$RNA@scale.data[markers,]))
for(g in markers){
  infoData[[g]]<-ave(infoData[[g]], infoData$cluster, FUN = mean)
}
meanData<-infoData[!duplicated(infoData$CRHBP),]
mat<-t(meanData[,2:ncol(meanData)])
colnames(mat)<-as.character(unique(infoData$cluster))
mat<-mat[,levels(infoData$cluster)]
mat2<-mat

colnames(mat1)<-paste0("Young ", colnames(mat1))
colnames(mat2)<-paste0("Elderly ", colnames(mat2))

mat<-cbind(mat1, mat2)
df<-data.frame(Condition = c(rep("Young", ncol(mat1)), 
                             rep("Elderly", ncol(mat2))),
               row.names = colnames(mat))
df$Condition<-factor(df$Condition, levels = c("Young", "Elderly"))
annot.col<-list(Condition=col.condition[1:2])
col<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(20)

pheatmap::pheatmap(mat, scale = "row", color = col, 
                   cluster_rows = F, cluster_cols = F, 
                   annotation_col = df, annotation_names_col = F,
                   annotation_colors = annot.col, width = 14, height = 9,
                   filename = "/home/mainciburu/scRNA/pics/figure1/heatmap_markers.png")
