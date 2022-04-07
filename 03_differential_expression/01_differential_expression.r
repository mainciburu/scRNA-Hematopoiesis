# Differential expression analysis for GSEA 
# Between celltype to characterize clusters
# between condition 
######################

library(Seurat)
library(foreach)
library(doParallel)

####### Input data #####################
## Repeat for 9 comparisons:
# young vs elderly / MDS1 / MDS2 / MDS3 / MDS4
# elderly vs MDS1 / MDS2 / MDS3 / MDS4

file1<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"
file2<-"/home/mainciburu/scRNA/MDS_paper/seurat_mds1.rds"
name1<-"senior"
name2<-"MDS1"
res.path<-"/home/mainciburu/scRNA/btwn_condition/"
btwn.celltype<-FALSE
btwn.condition<-TRUE
ncores<-4

seurat1<-readRDS(file1)
Idents(seurat1)<-seurat1$prediction
DefaultAssay(seurat1)<-"RNA"
seurat2<-readRDS(file2)
Idents(seurat2)<-seurat2$prediction
DefaultAssay(seurat2)<-"RNA"


########### Between celltype ################
if(btwn.celltype==TRUE){
	deg.seurat1<-FindAllMarkers(object = seurat1, test.use = "MAST", 
	                            assay = "RNA", only.pos=TRUE, logfc.threshold=0.1)
	res.file<-paste0(res.path, "deg_mast_celltype_", name1, ".rds")
	saveRDS(deg.seurat1, file = res.file)
	deg.seurat2<-FindAllMarkers(object = seurat2, test.use = "MAST", 
	                            assay = "RNA", only.pos=TRUE, logfc.threshold=0.1)
	res.file<-paste0(res.path, "deg_mast_celltype_", name2, ".rds")
	saveRDS(deg.seurat2, file = res.file)
}

########### Between condition ################

# Test for clusters where there are more than 25 cells in both conditions
if(btwn.condition==TRUE){	
	# Use only genes common in every individual
	g<-read.table("/home/mainciburu/scRNA/common_genes.txt", stringsAsFactors=F)$V1
	g<-g[g%in%rownames(seurat2)]
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
	  	                     assay = "RNA", test.use = "MAST", latent.vars = "Patient",
	  	                     logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
	    deg.tmp$gene<-rownames(deg.tmp)
	    deg.tmp$CellType<-cell
		}else{
		  deg.tmp<-NULL}
		print(dim(deg.tmp))
		return(deg.tmp)
	}
	res.file<-paste0(res.path, "deg_MAST_latent_var_", name1, "_", name2, ".Rdata")
	save(deg.condition, file = res.file)
}


## Expression proportion per sample
#res.prp<-list()
#for(ix in 1:length(celltypes)){
#	cell<-celltypes[ix]
#	print(paste0(cell, "-------------"))
#	dat.merge<-merge(x = subset(seurat1, idents = cell, features = g),     # seurat subset
#	                 y = subset(seurat2, idents = cell, features = g))
#	Idents(dat.merge)<-dat.merge$Condition
#	if(table(dat.merge$Condition)[1]>25 & table(dat.merge$Condition)[2]>25){
#		x<-dat.merge@assays$RNA@data
#		x[x>0]<-1
#		rs<-rowsum(t(as.matrix(x)), group = dat.merge$Patient)
#		prp<-apply(rs, 2, function(xx){xx/table(dat.merge$Patient)})
#		res.prp[[ix]]<-prp
#		names(res.prp)[ix]<-cell
#	}
#}

#tmp<-matrix(NA, ncol = length(g), nrow = 2)
#res.prp$ProB<-rbind(tmp, res.prp$ProB)
#rownames(res.prp$ProB)<-rownames(res.prp$HSC)

#load("scRNA/btwn_condition/deg_MAST_latent_var_young_senior.Rdata")

#for(ix in 1:length(deg.condition)){
#	if(nrow(res.prp[[ix]])==8){
#		i<-match(deg.condition[[ix]]$gene, colnames(res.prp[[ix]]))
#		deg.condition[[ix]]<-cbind(deg.condition[[ix]], t(res.prp[[ix]])[i,])
#		deg.condition[[ix]]<-deg.condition[[ix]][,c(1:7,11:15,8:10)]
#		colnames(deg.condition[[ix]])[8:15]<-paste0("pct.",c("Young1", "Young2", "Young3", "Young4", "Young5", "Elderly1", "Elderly2", "Elderly3"))
#	}
#}

#x<-do.call("rbind", deg.condition)
#write.csv(x, file = "scRNA/btwn_condition/deg_MAST_latent_var_young_senior.csv")

