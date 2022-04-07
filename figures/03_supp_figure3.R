######### supp figure 3 #################
# VlnPlots differentially expressed genes
#########################################

library(Seurat)
library(RColorBrewer)
library(ggsignif)
library(ggplot2)
source("/home/mainciburu/scRNA/colors.r")

################
young<-readRDS("scRNA/young/seurat_young_v3.rds")
senior<-readRDS("scRNA/senior/seurat_senior_v4.rds")

load("/home/mainciburu/scRNA/btwn_condition/deg_MAST_young_senior.Rdata")
deg<-do.call("rbind", deg.condition)

gg<-c("DDIT4", "IDS", "NFKBIA", "JUNB", "FOS", "FOSB", "ID2", "CD69",
      "MYC", "HMGN2", "CDK4", "UQCR11", "NDUFA13", "COX7A2", "PGK1", "TPI1")

x<-CreateSeuratObject(counts=young@assays$RNA@data[gg,], assay = "RNA",
                          min.cells = 0, min.features = 0, names.field = 1,
                          names.delim = "_", meta.data = young@meta.data)
x@assays$RNA@data<-x@assays$RNA@counts
x$Identity<-x$CellType2
y<-CreateSeuratObject(counts=senior@assays$RNA@data[gg,], assay = "RNA",
                      min.cells = 0, min.features = 0, names.field = 1,
                      names.delim = "_", meta.data = senior@meta.data)
y@assays$RNA@data<-y@assays$RNA@counts
Idents(y)<-"prediction"
y<-subset(y, idents = levels(young$CellType2))
y$prediction<-droplevels(y$prediction)
y$Identity<-y$prediction

dat<-merge(x = x, y = y)

Idents(dat)<-"Condition"
dat$Condition<-factor(dat$Condition, levels = c("Young", "Elderly"))
dat$Identity<-factor(dat$Identity, levels = levels(young$CellType2))
for(g in gg){
  df.sig<-deg[deg$gene==g & deg$p_val_adj<0.1,c(5,6,7)]
  df<-data.frame(Identity=levels(young$CellType2), 
                 pval=1,
                 significance=factor("NS", levels = c("NS", "*","**")))
  i<-match(df$Identity, df.sig$CellType)
  i<-i[!is.na(i)]
  df.sig<-df.sig[i,]
  df$pval[df$Identity%in%df.sig$CellType]<-df.sig$p_val_adj
  df$significance[df$pval<0.1]<-"*"
  df$significance[df$pval<0.01]<-"**"
  pdf(paste0("scRNA/figures/supp_figure3/vlnplot_", g, ".pdf"), 
      useDingbats = F, width = 15, height = 8)
  pp<-VlnPlot(dat, features = g, group.by = "Identity", split.by = "Condition", assay = "RNA", pt.size = 0) + 
    scale_fill_manual(values = c(Young="#D73027", Elderly ="#4575B4")) + ggtitle(label = g)
  n<-nrow(df)
  #pp<-pp + coord_cartesian(ylim = c(0, 6.5))
  ypos<-3
  pp<- pp + geom_signif(annotations = df$significance,
                        y_position = rep(ypos, n),
                        xmin = 1:n-0.3, xmax = 1:n + 0.3, tip_length = 0.01, 
                        size = 1, textsize = 8)
  pp<-pp + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(axis.text = element_text(size = 30), axis.title = element_text(size = 32), plot.title = element_text(size = 32)) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 32)) +
    theme(text = element_text(face = "plain")) +
    theme(axis.line.x =  element_line(size = 1.5), axis.ticks.x = element_line(size = 1.5)) + 
    theme(axis.line.y =  element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5))
  print(pp)
  dev.off()
}




