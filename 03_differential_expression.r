# Differential expression analysis for GSEA 
# Between celltype to characterize young clusters
# between condition to look for differences in comparable populations, using the glmnet predicted labels
  # only if there are more than 25 cells in both groups
######################

library(Seurat)
library(foreach)
library(doParallel)

####### Input data #####################
file1<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"
seurat1<-readRDS(file1)
file2<-"/home/mainciburu/scRNA/MDS/seurat_mds3.rds"
seurat2<-readRDS(file2)
Idents(seurat1)<-seurat1$prediction
DefaultAssay(seurat1)<-"RNA"
Idents(seurat2)<-seurat2$prediction
DefaultAssay(seurat2)<-"RNA"
btwn.celltype<-FALSE
btwn.condition<-TRUE


###### Between celltype #######################
if(btwn.celltype==TRUE){
	deg.seurat1<-FindAllMarkers(object = seurat1, test.use = "MAST", 
	                            assay = "RNA", only.pos=TRUE, logfc.threshold=0.1)
	res.file<-"/home/mainciburu/scRNA/btwn_celltype/deg_mast_celltype_senior.rds"
	saveRDS(deg.seurat1, file = res.file)
	deg.seurat2<-FindAllMarkers(object = seurat2, test.use = "MAST", 
	                            assay = "RNA", only.pos=TRUE, logfc.threshold=0.1)
	res.file<-"/home/mainciburu/scRNA/btwn_celltype/deg_mast_celltype_mds5.rds"
	saveRDS(deg.seurat2, file = res.file)
}

###### Between condition ######################

# Test for clusters where there are more than 25 cells in both conditions
if(btwn.condition==TRUE){	
	# Use only genes common in every individual
	g<-read.table("/home/mainciburu/scRNA/common_genes.txt", stringsAsFactors=F)$V1
	g<-g[g%in%rownames(seurat2)]
	ncores<-6
	registerDoParallel(cores=ncores)

	celltypes<-levels(Idents(seurat1))
	celltypes<-celltypes[celltypes%in%unique(seurat2$prediction)]
	celltypes<-celltypes[!celltypes%in%"not assigned"]

	ident.1<-unique(seurat1$Condition)
	ident.2<-unique(seurat2$Condition)

	deg.condition<-foreach(i = 1:length(celltypes), .packages = c("Seurat"))%dopar%{
		cell<-celltypes[i]
		print(paste0(cell, "-------------"))
		dat.merge<-merge(x = subset(seurat1, idents = cell, features = g),     # seurat subset
		                 y = subset(seurat2, idents = cell, features = g))
		Idents(dat.merge)<-dat.merge$Condition
		if(table(dat.merge$Condition)[1]>25 & table(dat.merge$Condition)[2]>25){
	  	deg.tmp<-FindMarkers(object = dat.merge, ident.1 = ident.1, ident.2 = ident.2, 
	  	                     assay = "RNA", test.use = "MAST", 
	  	                     logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
	    deg.tmp$gene<-rownames(deg.tmp)
	    deg.tmp$CellType<-cell
		}else{
		  deg.tmp<-NULL}
		print(dim(deg.tmp))
		return(deg.tmp)
	}

	res.file<-"/home/mainciburu/scRNA/btwn_condition/deg_MAST_senior_mds3.Rdata"
	save(deg.condition, file = res.file)
}

