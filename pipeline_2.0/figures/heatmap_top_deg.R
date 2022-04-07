# Heatmap top deg-------------------
library(ggplot2)
library(Seurat)
library(pheatmap)
library(RColorBrewer)

young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

load("/home/mainciburu/scRNA/btwn_condition/deg_MAST_young_senior_v4.Rdata")
i<-sapply(deg.condition, is.data.frame)
deg.condition<-deg.condition[i]
names(deg.condition)<-sapply(deg.condition, function(x){unique(x$CellType)})
pdf("/home/mainciburu/scRNA/figures_dec20/supp_figure3/heatmap_top_deg_MAST.pdf")
Idents(young)<-"CellType2"
Idents(senior)<-"prediction"
x<-intersect(young@assays$integrated@var.features, senior@assays$integrated@var.features)

# select celltypes with at least 10 cells per senior patient
t<-table(senior$prediction, senior$Patient)
rs<-rowSums(t>10)==3
rs<-names(rs)[rs==T]
deg.condition<-deg.condition[names(deg.condition)%in%rs]

for(i in 1:length(deg.condition)){
  celltype<-unique(deg.condition[[i]]$CellType)
  deg<-deg.condition[[i]]
  deg$pct.diff<-deg$pct.1 - deg$pct.2
  deg<-subset(deg, deg$p_val_adj<=1e-2)      # select significant genes
  deg<-deg[deg$gene%in%x,]        # select only variable genes
  g<-deg[order(deg$avg_logFC, decreasing = T),]       # order by FC
  rb<-c(grep("RPS", rownames(g)), grep("RPL", rownames(g)))  # remove ribosomal genes
  mt<-grep("MT-", rownames(g))                               # remove mitochondrial genes
  grm<-which(rownames(g)%in%c("AC009501.4", "CTD-2090I13.1", # remove sexual chromosome genes
                              "XIST", "PSMA2", "EIF1AY"))
  grm<-c(grm, rb, mt)
  g_up<-head(rownames(g)[-grm], 15)
  g_down<-tail(rownames(g)[-grm], 15)
  g_top<-c(g_up, g_down)
  g.merge<-merge(x = subset(young, idents = celltype),
                 y = subset(senior, idents = celltype))
  g.merge<-ScaleData(g.merge, assay = "RNA")
  dat<-data.frame(Patient=g.merge$Patient,
                  row.names = rownames(g.merge@meta.data))
  dat<-cbind(dat, t(g.merge@assays$RNA@scale.data[g_top,]))
  dat$Patient<-plyr::mapvalues(x = dat$Patient, from = unique(dat$Patient),
                               to = c("Young1", "Young2", "Young3", "Young4","Young5",
                                      "Elderly1", "Elderly2", "Elderly3"))
  dat$Patient<-factor(dat$Patient, 
                      levels = c("Young1", "Young2", "Young3", "Young4", "Young5",
                                 "Elderly1", "Elderly2", "Elderly3"))
  for(gg in g_top){
    dat[[gg]]<-ave(dat[[gg]], dat$Patient, FUN = mean)
  }
  dd<-apply(dat[,2:30], 2, function(x){sum(!duplicated(x))})
  jj<-as.numeric(which(dd==nlevels(dat$Patient))[1]) + 1
  meanData<-dat[!duplicated(dat[,jj]),]
  mat<-t(meanData[,2:ncol(meanData)])
  colnames(mat)<-as.character(unique(dat$Patient))
  col<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(20)
  col<-viridis::viridis(n = 30)
  df<-data.frame(Condition = c(rep("Young", 5), rep("Elderly", 3)), 
                 row.names = colnames(mat))
  annot.col<-list(Condition=col.condition[1:2])
  print(pheatmap::pheatmap(mat, scale = "row", color = col, 
                           cluster_rows = F, cluster_cols = F, 
                           annotation_col = df, annotation_names_col = F, 
                           annotation_legend = F, treeheight_row = 0, 
                           fontsize = 12, cellwidth = 30, main = celltype, annotation_colors = annot.col))
}
dev.off()
