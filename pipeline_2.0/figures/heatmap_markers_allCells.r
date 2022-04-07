young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v3_with_predictions.rds")

# Heatmaps ---------------------------- 
hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
cc<-c("CDC20", "TOP2A")
gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
granul<-c("ELANE", "AZU1", "CEBPA", "CEBPE", "CST7")
mono<-c("LYZ", "CSTA", "CD14")
dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4", "IGKC", "SPIB", "LILRA4", "GZMB", "BLINK")
t<-c("JCHAIN", "PRSS2", "TXNP", "LTB")
nk<-c("TSC22D1", "CXXC5", "HOXA9", "HOXA10")
clp<-c("IL7R", "DNTT")
prob<-c("VPREB1", "EBF1", "CD79A", "CD79B", "IGHM")
mep<-c("NFE2", "HFS1", "TAL1")
mk<-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")
GX<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, prob, t, nk, mep, mk, ery, baso)

# Young celltype ----------------------------------------------------------
Counts<-young@assays$RNA@counts[,rownames(young@meta.data)]
Counts<-as.matrix(Counts)
Counts<-(t(t(Counts)/rowSums(t(Counts))))*1e6
Data<-log10(Counts+1)

colannotations<-data.frame(CellType=young@meta.data$CellType2)
rownames(colannotations)<-colnames(Data)
colannotations<-colannotations[order(colannotations$CellType),,drop=F]

iix<-match(GX,rownames(Data))
GXi<-GX[-which(is.na(iix))]
hmcols<-colorRampPalette(c("#ffffff","#fcf1de","#FFC75F","#FF9671","#FF6F91","#845EC2"))

pheatmap(Data[GXi,rownames(colannotations)],annotation_col=colannotations,show_colnames = FALSE,show_rownames = TRUE,cluster_cols = FALSE,
         cluster_rows = FALSE,color=hmcols(10),annotation_colors = list(CellType=col.young.v3),
         border_color = "black",legend=T,annotation_legend = F,annotation_names_row = F,annotation_names_col = F,
         filename="/home/mainciburu/scRNA/figures_nov20/figure1/heatmap_markers_young.pdf", width=8, height=12)

# senior predictions ----------------------------------------------------
Counts<-senior@assays$RNA@counts[,rownames(senior@meta.data)]
Counts<-as.matrix(Counts)
Counts<-(t(t(Counts)/rowSums(t(Counts))))*1e6
Data<-log10(Counts+1)

colannotations<-data.frame(CellType=senior@meta.data$prediction)
rownames(colannotations)<-colnames(Data)
colannotations<-colannotations[order(colannotations$CellType),,drop=F]

iix<-match(GX,rownames(Data))
GXi<-GX[-which(is.na(iix))]
hmcols<-colorRampPalette(c("#ffffff","#fcf1de","#FFC75F","#FF9671","#FF6F91","#845EC2"))

col<-c(col.young.v3, "not assigned"="grey")

pheatmap(Data[GXi,rownames(colannotations)],annotation_col=colannotations,show_colnames = FALSE,show_rownames = TRUE,cluster_cols = FALSE,
         cluster_rows = FALSE,color=hmcols(10),annotation_colors = list(CellType=col),
         border_color = "black",legend=T,annotation_legend = F,annotation_names_row = F,annotation_names_col = F,
         filename="/home/mainciburu/scRNA/figures_nov20/figure1/heatmap_markers_senior.pdf", width=8, height=12)
