library(Seurat)
library(ggplot2)

young<-readRDS("/home//mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v3.rds")

hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
cc<-c("CDC20", "TOP2A")
gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
granul<-c("ELANE", "AZU1", "CEBPA", "CEBPE", "CST7")
mono<-c("LYZ", "CSTA", "CD14")
dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4", "IGKC", "SPIB", "LILRA4", "GZMB", "BLINK")
t<-c("JCHAIN", "PRSS2", "TXNP", "LTB")
nk<-c("CXXC5", "HOXA9", "HOXA10")
clp<-c("IL7R", "DNTT")
prob<-c("VPREB1", "EBF1", "CD79A", "CD79B", "IGHM")
mep<-c("NFE2", "HFS1", "TAL1")
mk<-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")

GX<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, prob, t, nk, mep, mk, ery, baso)

GX<-data.frame(Gene=GX,                                                   # Normalized counts for marker genes
               Senior=match(GX,rownames(senior@assays$RNA@data)),
               Young=match(GX,rownames(young@assays$RNA@data)), stringsAsFactors = F)
GX<-GX[-which(is.na(GX$Senior) & is.na(GX$Young)),]

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

ggsave(file="/home/mainciburu/scRNA/figures_nov20/figure1/DotPlot_YoungSenior.pdf",
       height = 10)
