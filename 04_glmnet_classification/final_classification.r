library(knitr)
library(glmnet)
library(caret)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(gridExtra)



########### input data ########################
args=(commandArgs(TRUE))

#train.path<-args[1]
test.path<-args[1]            # seurat object for prediction
res.path<-args[2]             # folder to store results
celltypes.file<-args[3]       # text file with every cell types
orig.group<-args[4]           # original classification, stored as a column in metadata, for the test object
colors<-args[5]               # colors for the plot
reduction<-args[6]            # reduction used for visualization

print("Imported arguments---------")
print(paste0("test.path = ", test.path))
print(paste0("res.path = ", res.path))
print(paste0("celltype = ", celltype))
print(paste0("orig.group = ", orig.group))
print(paste0("colors = ", colors))
print(paste0("reduction = ", reduction))

source("/home/mainciburu/scRNA/colors.r")
colors<-col.pred

# Load seurat
print("Loading data-------------")
#seurat.train<-readRDS(train.path)
seurat.test<-readRDS(test.path)

###########################################################################################

# celltype names
celltypes<-read.table(file = file.path(celltypes.file), header = F, col.names = "celltypes", sep = "\n", stringsAsFactors = F)

# load predictions from every binary model
# store the average probability (and the total number of positives)
pred.avg.prob<-data.frame(row.names = colnames(seurat.test))
pred.pos.sum<-data.frame(row.names = colnames(seurat.test))

for(cell in celltypes$celltypes){
  pred.tmp<-read.table(file = file.path(paste0(res.path, cell, "_predictions.txt")))
  pred.avg.prob[,cell]<-pred.tmp$avg.prob
  pred.pos.sum[,cell]<-pred.tmp$pos.sum
}

# extract maximun probability and its label
pred.avg.prob$max<-apply(pred.avg.prob[,1:ncol(pred.avg.prob)], 1, max)
pred.avg.prob$max.label<-apply(pred.avg.prob[,1:(ncol(pred.avg.prob)-1)], 1,
							   function(x){
							   		return(colnames(pred.avg.prob)[which(x==max(x))])
							   		}
							   	)

# if max avg.prob <= 0.5 => assign no label
pred.avg.prob$max.label[pred.avg.prob$max<=0.5]<-"not assigned"

#pred.pos.sum$max<-apply(pred.pos.sum[,1:ncol(pred.pos.sum)], 1, max)
#pred.pos.sum$max.label<-apply(pred.pos.sum[,1:(ncol(pred.pos.sum)-1)], 1, function(x){
#  m<-colnames(pred.pos.sum)[which(x==max(x))]
#  if(length(m)==1){return(m)}
#  else{return("tie")}
#})

# if max pos.sum <= 5 --> assign no label
#pred.pos.sum$max.label[pred.pos.sum$max<=5]<-"not assigned"

# Store final prediction in the Seurat objects
seurat.test$prediction<-pred.avg.prob$max.label
seurat.test$prediction.prob<-pred.avg.prob$max


seurat.test$prediction<-factor(seurat.test$prediction, 
                                     levels = c(CellTypes, "not assigned"))

# Save the object with predictions
new.path<-gsub(pattern = ".rds", replacement = "_with_predictions.rds", x = test.path)
saveRDS(seurat.test, file = new.path)

# generate pdf plots
pdf(paste0(res.path, "final_predictions.pdf"), width = 15, height = 10, useDingbats = F)
# plot final predictions
DimPlot(seurat.test, reduction = reduction, group.by = "prediction", cols = colors) + ggtitle("Final prediction")
FeaturePlot(seurat.test, features = "prediction.prob", reduction = reduction, 
            cols = rev(brewer.pal(11, "RdYlBu"))) + ggtitle("Probability")

# table original groups vs glmnet final predictions
t<-table(seurat.test@meta.data[,orig.group], seurat.test$prediction)
frame()
gridExtra::grid.table(t)
dev.off()

