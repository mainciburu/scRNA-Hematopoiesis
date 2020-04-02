# Create model for a specific cell type 
# Cell type is specified in the .sbs file
# output: matrix with predictions per cell and iteration + average prediction
#         pdf plots of prediction probabilities + binary labels + table original classification vs predictions

library(knitr)
library(glmnet)
library(ROCR)
library(pROC)
library(caret)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(gridExtra)


########## Functions ################################


#####################################################



########## Input data ###############################
args=(commandArgs(TRUE))

train.path<-args[1]     # seurat object for training and validation
test.path<-args[2]      # seurat object for prediction
res.path<-args[3]       # folder to store results
celltype<-args[4]       # specific cell type for the model
orig.group<-args[5]     # original classification, stored as a column in metadata, for the test object
reduction<-args[6]      # reduction used for visualization

print("Imported arguments---------")
print(paste0("train.path = ", train.path))
print(paste0("test.path = ", test.path))
print(paste0("res.path = ", res.path))
print(paste0("celltype = ", celltype))
print(paste0("orig.group = ", orig.group))
print(paste0("reduction = ", reduction))

print("Loading data-------------")
seurat.train<-readRDS(train.path)
seurat.test<-readRDS(test.path)

#####################################################

# Cycling LMPP => Cycling_LMPP
levels(seurat.train$CellType)<-c("HSC", "LMPP", "Cycling_LMPP", levels(seurat.train$CellType)[4:12])
Idents(seurat.train)<-"CellType"

# select features: variable genes present in both datasets
x<-seurat.train@assays$RNA@var.features
x<-x[x%in%rownames(seurat.test)]

# Iteration number
r<-10

print(paste0("###Build ", celltype, " model"))

# Create meta.data column with binary cellType label
seurat.train@meta.data[,celltype]<-0
seurat.train@meta.data[,celltype][seurat.train$CellType == celltype]<-1

# Empty matrix to store future predictions 
pred.tmp<-matrix(nrow = nrow(seurat.test@meta.data), ncol = r)
rownames(pred.tmp)<-rownames(seurat.test@meta.data)
colnames(pred.tmp)<-paste0("Round_", 1:r)


############ START LOOP ################################
for(z in 1:r){
  print(paste0("Round ", z, "----------"))
  # select cells: same number of celltype and non_celltype
  y<-seurat.train@meta.data[,celltype][seurat.train@meta.data[,celltype]==1]
  names(y)<-rownames(seurat.train@meta.data)[seurat.train@meta.data[,celltype]==1]
  y2<-seurat.train@meta.data[,celltype][seurat.train@meta.data[,celltype]==0]
  names(y2)<-rownames(seurat.train@meta.data)[seurat.train@meta.data[,celltype]==0]
  i<-sample(x = 1:length(y2), size = length(y), replace = F)
  y2<-y2[i]
  y<-c(y, y2)

  # Divide y: 75% train, 25% validation
  keep_train<-sample(x = 1:length(y), size = 0.75*length(y), replace = F)
  y_train<-y[keep_train]
  y_validation<-y[-keep_train]
  
  # Build train and validation matrix (cells x genes)
  x_train<-t(as.matrix(seurat.train@assays$RNA@data[x, names(y_train)]))
  x_validation<-t(as.matrix(seurat.train@assays$RNA@data[x, names(y_validation)]))
  
  # Fit model: try alpha =0.25, 0.5, 0.75, 1
  # Same cell foldid for all the models
  
  foldid<-sample(1:10, size = length(y_train), replace = T)
  print("    fit cv1")
  cv1<-cv.glmnet(x = x_train, y = y_train, foldid = foldid, alpha = 1, family = "binomial", type.measure = "auc")
  print("    fit cv.75")
  cv.75<-cv.glmnet(x = x_train, y = y_train, foldid = foldid, alpha = 0.75, family = "binomial", type.measure = "auc")
  print("    fit cv.5")
  cv.5<-cv.glmnet(x = x_train, y = y_train, foldid = foldid, alpha = 0.5, family = "binomial", type.measure = "auc")
  print("    fit cv.25")
  cv.25<-cv.glmnet(x = x_train, y = y_train, foldid = foldid, alpha = 0.25, family = "binomial", type.measure = "auc")
  print("    fit cv.1")
  cv.1<-cv.glmnet(x = x_train, y = y_train, foldid = foldid, alpha = 0.1, family = "binomial", type.measure = "auc")
  
  
  # Test all models on the validation set
  print("    Validate")
  cv<-list(cv1=cv1, cv.75=cv.75, cv.5=cv.5, cv.25=cv.25, cv.1=cv.1)
  eval<-data.frame(TPR=NA, TNR=NA, FPR=NA, FNR=NA, AUC=NA, Vars=NA)
  for(i in 1:length(cv)){
    pred_val<-predict(cv[[i]], x_validation, s="lambda.1se", type="response")     # probabilities to calculate auc
    auc_val<-pROC::roc(y_validation, pred_val[,1])
    auc_val<-as.numeric(pROC::auc(auc_val))
    eval[i,5]<-auc_val 
    pred_val<-predict(cv[[i]], x_validation, s="lambda.1se", type="class")        # classes
    pred_val<-as.factor(pred_val[,1])
    t<-caret::confusionMatrix(data = pred_val, as.factor(y_validation), positive = "1")$table
    eval[i,1]<-t[2,2]/(t[2,2]+t[2,1])    # TP
    eval[i,2]<-t[1,1]/(t[1,1]+t[1,2])    # TN
    eval[i,3]<-t[2,1]/(t[2,2]+t[2,1])    # FP
    eval[i,4]<-t[1,2]/(t[2,2]+t[2,1])    # FN
    v<-coef(cv[[i]], s = cv[[i]]$lambda.1se)    # number of selected variables
    v<-v[,1]
    eval[i,6]<-sum(v!=0)
  }
  
  rownames(eval)<-c("cv1","cv.75","cv.5","cv.25","cv.1")
  
  # Choose best alpha

  # 1) vars in [20,150]
    # ** there might be no model with > 20 vars or < 150 vars
  if(sum(eval$Vars>=20)>=1&sum(eval$Vars<=150)>=1){
    best<-eval[eval$Vars>=20&eval$Vars<=150,]
    best<-round(best, 3)
  }else{
  	best<-eval
 	  best<-round(best, 3)
  }
  # 2) max auc
  if(nrow(best)>1){
    i<-max(best$AUC)
    best<-best[best$AUC==i,]
  }
  # 3) min FPR
  if(nrow(best)>1){
    i<-min(best$FPR)
    best<-best[best$FPR==i,]
  }
  # 4) min FNR
  if(nrow(best)>1){
    i<-min(best$FNR)
    best<-best[best$FNR==i,]
  }
  # 5) max variables
  if(nrow(best)>1){
    i<-max(best$Vars)
    best<-best[best$Vars==i,]
  }
  # if still more than 1 model => choose the first one
  if(nrow(best)>1){
    best<-best[1,]
  }
  best<-rownames(best)
  print(paste0("    Optimal alpha: ", best))
  best<-which(names(cv)==best)
  glmnet.opt<-cv[[best]]
  
  # Test best model on the test seurat object
  print("    Test")
  test<-TRUE
  
  if(test){
    x_test<-t(as.matrix(seurat.test@assays$RNA@data[x,]))
    pred_test<-predict(glmnet.opt, x_test, type="response")
    pred.tmp[,z]<-pred_test[,1]     # store probabilities for one iteration
  }
}  
print("Loop ended---------")
####################### END LOOP #################################################

print("Summarize predictions-------------")
pred.tmp<-as.data.frame(pred.tmp)
pred.tmp$avg.prob<-rowMeans(pred.tmp)    # average probabilities among all iterations
pred.tmp$sd.prob<-apply(pred.tmp[,1:r], 1, sd)
pred.tmp$pos.sum<-rowSums(pred.tmp[,1:r]>0.5)     # total number of positives among all iterations
pred.tmp$final.res<-0
pred.tmp$final.res[pred.tmp$avg.prob>0.5]<-1

# Save results
write.table(x = pred.tmp, file = paste0(res.path, celltype, "_predictions.txt"))

print("Generate plots--------")
# plot avg probability and final result on test object
seurat.test$pred<-pred.tmp$avg.prob
p1<-FeaturePlot(seurat.test, features = "pred", reduction = reduction, 
                cols = rev(brewer.pal(11, "RdYlBu"))) + ggtitle("Probability")
seurat.test$pred<-pred.tmp$final.res
p2<-DimPlot(seurat.test, group.by =  "pred", reduction = reduction, 
            cols = c("grey", "red")) + ggtitle("Class")
pp<-plot_grid(p1, p2, ncol = 2)

# Glmnet results vs original groups
eval<-table(seurat.test@meta.data[,orig.group], pred.tmp$final.res)
colnames(eval)<- c("Negative", "Positive")
tt<-eval


# Export pdf
pdf(paste0(res.path, celltype, "_results.pdf"), width = 15, height = 10, useDingbats = FALSE)
pp
frame()
gridExtra::grid.table(tt)
dev.off()
