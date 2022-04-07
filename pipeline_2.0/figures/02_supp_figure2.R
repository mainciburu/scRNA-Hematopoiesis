########## supplementary figure 2 ################
# UMAP every prediction score
# UMAP original greenleaf groups
# UMAP predicted labels
# AUC
##################################################

library(knitr)
library(glmnet)
library(caret)
library(AUC)
library(ROCR)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

source("/home/mainciburu/scRNA/colors.r")

# UMAP prediction scores----------------------
celltypes<-read.table(file = "/home/mainciburu/scRNA/pipeline_2.0/04_glmnet_classification/celltypes.txt",
                      header = F, col.names = "celltypes", sep = "\n", stringsAsFactors = F)
cols<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(30)
pred.avg.prob<-data.frame(row.names = colnames(senior))
for(cell in celltypes$celltypes){
  pred.tmp<-read.table(file = file.path(paste0("/home/mainciburu/scRNA/senior/glmnet_v4/", cell, "_predictions.txt")))
  i<-match(rownames(senior@meta.data), rownames(pred.tmp))
  senior$pred<-pred.tmp$avg.prob[i]
  pdf(file = paste0("/home/mainciburu/scRNA/figures_mar21/supp_figure2/UMAP_elderly_", cell, "_score.pdf"), useDingbats = F,
      width = 6, height = 5)
  print(FeaturePlot(senior, reduction = "umap.int", features = "pred",
              pt.size = 0.5, cols = cols) + ggtitle(cell) +
    labs(colour = "Prediction score")  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2"))
  dev.off()
}

# UMAP Greenleaf --------------------------------
gr<-readRDS("/home/mainciburu/scRNA/greenleaf/seurat_greenleaf_cd34.rds")
# Original groups
pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure2/UMAP_greenleaf_original.pdf", useDingbats = F,
    width = 10, height = 10)
DimPlot(gr, group.by = "BioClassification", cols = col.greenleaf, pt.size = 0.5, label = F) + 
  labs(colour = "Original Groups")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

# Predicted labels
pdf(file = "/home/mainciburu/scRNA/figures_mar21/supp_figure2/UMAP_greenleaf_predictions.pdf", useDingbats = F,
    width = 10, height = 10)
col<-c(col.young.v3[names(col.young.v3)%in%levels(gr$prediction)], "not assigned"="grey")
DimPlot(gr, group.by = "prediction", cols = col, pt.size = 0.5, label = F) + 
  labs(colour = "Predicted Cell Type")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2") 
dev.off()

# AUC -----------------------------------
# List with corresponding labels
match.labels<-list(HSC="HSC", 
                   LMPP="CMP.LMPP",
                   GMP_Granulocytes = c("GMP", "GMP.Neut"),
                   GMP = c("GMP", "GMP.Neut"),
                   Monocytes = "CD14.Mono.2",
                   pDC = "pDC",
                   CLP = c("CLP.1", "CLP.2"),
                   ProB = "Pre.B",
                   T_NK = c("CLP.1", "CLP.2"),
                   MEP = "Early.Eryth",
                   Megakaryocytes="Early.Eryth",
                   Erythroid_early = c("Early.Eryth", "Late.Eryth"),
                   Erythroid_late = c("Early.Eryth", "Late.Eryth"),
                   Basophils = "Early.Baso")

# Load binary models results
pred.avg.prob<-data.frame(row.names = colnames(gr))

celltypes<-read.table(file = "/home/mainciburu/scRNA/pipeline_2.0/04_glmnet_classification/celltypes.txt",
                      header = F, col.names = "celltypes", sep = "\n", stringsAsFactors = F)
for(cell in celltypes$celltypes){
  pred.tmp<-read.table(file = file.path(paste0("/home/mainciburu/scRNA/greenleaf/glmnet_v3/", 
                                               cell, "_predictions.txt")))
  pred.avg.prob[,cell]<-pred.tmp$avg.prob
}

pred.avg.prob$max<-apply(pred.avg.prob[,1:ncol(pred.avg.prob)], 1, max)
pred.avg.prob$max.label<-apply(pred.avg.prob[,1:(ncol(pred.avg.prob)-1)], 1, function(x){return(colnames(pred.avg.prob)[which(x==max(x))])})

# if max avg.prob <= 0.5 --> assign no label
pred.avg.prob$max.label[pred.avg.prob$max<=0.5]<-"not assigned"

# auc for every model
res<-apply(celltypes, 1, function(cell){
  pred<-pred.avg.prob[[cell]]       # vector prediction score per cell
  names(pred)<-rownames(pred.avg.prob)
  ix<-match.labels[[cell]]
  y<-rep(0, length(pred))     # vector correct labels -> 0 = negative, 1 = positive
  names(y)<-names(pred)
  y[gr$BioClassification%in%ix]<-1
  y<-factor(y, levels = c(0, 1))
  pred<-prediction(predictions = pred, labels = y, label.ordering = c("0", "1"))
  perf<- performance(pred,"tpr","fpr")
  auc1<-performance(pred,"auc")@y.values[[1]]
  df<-data.frame(fp=perf@x.values[[1]], tp=perf@y.values[[1]], AUC = auc1, celltype=cell)
  return(df)
  #res[[cell]]<-df
  #ggplot(df, aes(fp, tp)) + geom_line()
  #plot(perf@x.values[[1]],perf@y.values[[1]],col="black",type="l",lwd=2,bty="n",
  #     xlab="False positive rate",ylab="True positive rate", main = cell, sub = paste0("AUC = ", round(AUC1, 3)))
  
}
)
res<-do.call(rbind, res)
AUC<-res[!duplicated(res$AUC),]
AUC$fp<-0.75
AUC$tp<-0.25
pdf("/home/mainciburu/scRNA/figures_mar21/supp_figure2/roc_greenleaf.pdf", width = 12)
ggplot(res, aes(fp, tp)) + geom_line() + facet_wrap(~ celltype, ncol = 3) + 
  geom_text(data = AUC, label = paste0("AUC = ", round(AUC$AUC, 3))) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold"))
dev.off()
