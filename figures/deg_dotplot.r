### DEG dotplot
# color by mean expression
# size by pct of cells expressing
# divided by condition
library(ggplot2)
library(RColorBrewer)
library(Seurat)

load("/home/mainciburu/scRNA/btwn_condition/deg_MAST_young_senior.Rdata")
i<-sapply(deg.condition, is.data.frame)
deg.condition<-deg.condition[i]
names(deg.condition)<-sapply(deg.condition, function(x){unique(x$CellType)})

g_top<-c()
for(i in 1:length(deg.condition)){
  celltype<-unique(deg.condition[[i]]$CellType)
  deg<-deg.condition[[i]]
  deg$pct.diff<-deg$pct.1 - deg$pct.2
  deg<-subset(deg, deg$p_val_adj<=1e-2)
  g<-deg[order(deg$avg_logFC, decreasing = T),]
  rb<-c(grep("RPS", rownames(g)), grep("RPL", rownames(g)))  # remove ribosomal genes
  mt<-grep("MT-", rownames(g))  # remove mitocondrial genes
  grm<-which(rownames(g)%in%c("AC009501.4", "CTD-2090I13.1", # remove sexual chromosome genes
                              "XIST", "PSMA2", "EIF1AY"))
  grm<-c(grm, rb, mt)
  g_up<-head(rownames(g)[-grm], 5)
  g_down<-tail(rownames(g)[-grm], 5)
  g_top<-unique(c(g_top, g_up, g_down))
}

young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

cells<-names(deg.condition)
Idents(young)<-"CellType2"
Idents(senior)<-"prediction"
young<-subset(young, idents=cells)
senior<-subset(senior, idents=cells)

# A) Seurat dotplot
DefaultAssay(young)<-"RNA"
DefaultAssay(senior)<-"RNA"
P1<-DotPlot(young, assay = "RNA",features = g_top,dot.min = 0.2,group.by = "CellType2",
            col.min = -2,col.max = 2,col=c("#f7e69c","red"))

P2<-DotPlot(senior,assay = "RNA",features = g_top,dot.min = 0.2,group.by = "prediction",
            col.min = -2,col.max = 2,col=c("#f7e69c","red"))

YoungData<-P1$data
YoungData$Condition<-"Young"

SeniorData<-P2$data
SeniorData$Condition<-"Senior"

PlotData<-rbind(YoungData,SeniorData)
PlotData$Condition<-factor(PlotData$Condition,levels = c("Young","Senior"))
PlotData$features.plot<-factor(PlotData$features.plot,levels = rev(levels(PlotData$features.plot)))
cols<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(30)
ggplot(PlotData,aes(x=features.plot,y=Condition,colour=avg.exp.scaled)) +
  geom_point(aes(size=pct.exp))+facet_grid(~id)+
  scale_colour_gradientn(colours = cols)+theme_classic()+coord_flip()

#low = "#f7e69c",high = "red"
#low = "#313695",high = "#A50026"

# B) ggplot2 dotplot
df<-data.frame()
for(i in 1:length(cells)){
  a<-young@assays$RNA@data[g_top, young$CellType2==cells[i]]
  a<-as.matrix(a)
  a<-t(scale(t(a)))
  a<-rowMeans(a)
  b<-senior@assays$RNA@data[g_top, senior$prediction==cells[i]]
  b<-as.matrix(b)
  b<-t(scale(t(b)))
  b<-rowMeans(b)
  ix<-na.omit(match(g_top,deg.condition[[i]]$gene))
  pct1<-deg.condition[[i]]$pct.1[ix]
  pct2<-deg.condition[[i]]$pct.2[ix]
  a<-data.frame(Gene=g_top, CellType=cells[i], Condition="Young", MeanExpression=a, 
                pctExpression=pct1)
  b<-data.frame(Gene=g_top, CellType=cells[i], Condition="Elderly", MeanExpression=b, 
                pctExpression=pct2)
  df<-rbind(df, a, b)
}

cols<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(30)
pdf("dot.pdf", height = 15)
ggplot(df, aes(Condition, Gene)) + 
  geom_point(aes(fill = MeanExpression, shape=21, size = pctExpression)) + scale_shape_identity() +
  facet_grid(~CellType) +
  theme_bw() + scale_fill_gradientn(colours = cols)
dev.off()