# Extract from Seurat object
  # loom object
  # Normalized expression matrix
  # Integrated expression matrix
  # umap coordinates


library(loomR)
library(Seurat)
#library(hdf5r)

### young and senior
seurat.file<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"
seurat<-readRDS(seurat.file)

res.path<-"/home/mainciburu/scRNA/palantir/data/"
res.name<-"senior"
assay<-c("integrated", "RNA")

####### .loom objects #########
# Remove meta.data columns with NA
seurat@meta.data<-seurat@meta.data[,colSums(is.na(seurat@meta.data))==0]

# Create .loom object with integrated expression matrix (with every common genes)
if("integrated"%in%assay){
	#g<-seurat@assays$integrated@var.features
	dat<-CreateSeuratObject(counts = seurat@assays$integrated@data, assay = "RNA", meta.data = seurat@meta.data)
	dat@assays$RNA@data<-dat@assays$RNA@counts
	dat@assays$RNA@var.features<-rownames(dat@assays$RNA@counts)
	dat.loom<-as.loom(dat, filename = paste0(res.path, res.name, "_int.loom"))

	dat.loom$close_all()
}

# Create .loom object with original normalized expression matrix 
# include all common genes
if("RNA"%in%assay){
	g<-read.table("/home/mainciburu/scRNA/common_genes.txt")$V1
	dat<-CreateSeuratObject(counts = seurat@assays$RNA@data[g,], assay = "RNA", meta.data = seurat@meta.data)
	dat@assays$RNA@data<-dat@assays$RNA@counts
	dat@assays$RNA@var.features<-rownames(dat@assays$RNA@counts)
	dat.loom<-as.loom(dat, filename = paste0(res.path, res.name, "_norm.loom"))
	dat.loom$close_all()
}

####### umap ###########
# Get integrated umap coordinates
umap<-seurat@reductions$umap.int@cell.embeddings
write.table(x = umap, file = paste0(res.path, res.name, "_umap.txt"), sep = "\t", quote = F, col.names = F, row.names = T)


######## csv expression matrix #############
# integrated matrix with every common gene
if("integrated"%in%assay){
  x<-seurat@assays$integrated@data
  write.csv(x = x, file = paste0(res.path, res.name, "_int_mat.csv"))
}

# Full original normalized matrix
# include all common genes
g<-read.table("/home/mainciburu/scRNA/common_genes.txt")$V1
x<-seurat@assays$RNA@data[g,]
write.csv(x = x, file = paste0(res.path, res.name, "_norm_mat_full.csv"))


#################### MDS5 ####################
seurat.file<-"/home/mainciburu/scRNA/MDS/seurat_mds5.rds"
seurat<-readRDS(seurat.file)

res.path<-"/home/mainciburu/scRNA/palantir/data/"
res.name<-"mds5"
seurat@meta.data<-seurat@meta.data[,colSums(is.na(seurat@meta.data))==0]

## loom object
g<-read.table("/home/mainciburu/scRNA/common_genes.txt")$V1
g<-g[g%in%rownames(seurat)]
dat<-CreateSeuratObject(counts = seurat@assays$RNA@data[g,], assay = "RNA", meta.data = seurat@meta.data)
dat@assays$RNA@data<-dat@assays$RNA@counts
dat@assays$RNA@var.features<-rownames(dat@assays$RNA@counts)
dat.loom<-as.loom(dat, filename = paste0(res.path, res.name, "_norm.loom"))
dat.loom$close_all()

## umap coordinates
umap<-seurat@reductions$umap@cell.embeddings
write.table(x = umap, file = paste0(res.path, res.name, "_umap.txt"), sep = "\t", quote = F, col.names = F, row.names = T)

## csv
x<-seurat@assays$RNA@data[g,]
write.csv(x = x, file = paste0(res.path, res.name, "_norm_mat_full.csv"))



