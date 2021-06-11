#### Integrate samples from different conditions (Yong + Senior or Young + MDS)
# We only want the umap coordinates from this integration
# We add the new coordinates to the already created objects, integrated for each condition

library(Seurat)
library(umap)
library(dplyr)
library(RColorBrewer)
library(cluster)
library(future)

# Load integrated objects
seurat1.file<-"/home/mainciburu/scRNA/young/seurat_young_v3.rds"
seurat2.file<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"

seurat1<-readRDS(seurat1.file)
seurat2<-readRDS(seurat2.file)

# Merge the 2 objects. Use the original (RNA) data
DefaultAssay(seurat1)<-"RNA"
DefaultAssay(seurat2)<-"RNA"

AllData<-merge(seurat1, seurat2)

# split by patient
AllData<-SplitObject(object = AllData, split.by = "Patient")

# Normalize, scale, find variable genes (per patient)
#for(i in 1:5){
#  AllData[[i]]<-NormalizeData(AllData[[i]], scale.factor = 10000)
#  AllData[[i]]<-ScaleData(AllData[[i]], features = rownames(AllData[[i]]))
#  AllData[[i]]<-FindVariableFeatures(AllData[[i]], selection.method = "vst", nfeatures = 3000)
#}

# Activate parallelization
plan("multiprocess", workers = 8)
options(future.globals.maxSize= 20000000000)   ## change

# Integrate - 50 dims and 3000 features
#nfeatures<-3000
ndim<-50
#int.features<-SelectIntegrationFeatures(object.list = AllData, nfeatures = nfeatures)
anchors<-FindIntegrationAnchors(object.list = AllData, dims = 1:ndim)
# integrate common genes 
g<-read.table("/home/mainciburu/scRNA/common_genes.txt", header = F, stringsAsFactors = F)
g<-g$V1
seurat.int<-IntegrateData(anchors, dims = 1:ndim, features = g)

# Recalculate cell Cycle score 
DefaultAssay(seurat.int)<-"RNA"
seurat.int<-CellCycleScoring(seurat.int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
# scale integrated data and regress cc effect
DefaultAssay(seurat.int)<-"integrated"
seurat.int<-ScaleData(seurat.int, vars.to.regress = c("S.Score", "G2M.Score"))

# PCA
seurat.int<-RunPCA(seurat.int, npcs = ndim)
#DimPlot(seurat.int, dims = 1:2, reduction = "pca", group.by ="Patient", pt.size = 0.1)
# decide dimensionality
ElbowPlot(seurat.int, ndims = ndim)

# umap
npcs<-30
reduction.name<-"umap_integrated"
seurat.int<-RunUMAP(object = seurat.int, assay = "integrated", dims = 1:npcs, 
                    reduction.name = reduction.name, seed.use = 123)
pdf("umap_30pc.pdf")
DimPlot(seurat.int, reduction = "umap_integrated", dims = 1:2, group.by = "CellType2", pt.size = 0.1)
dev.off()

# Save umap coordinates as txt
umap<-seurat.int@reductions$umap_integrated@cell.embeddings
umap.file<-"/home/mainciburu/scRNA/umap_young_senior.txt"
write.table(umap, file = umap.file, col.names = F)

# Add coordinates to the original objects
umap1<-umap[match(colnames(seurat1), rownames(umap)),]
umap2<-umap[match(colnames(seurat2), rownames(umap)),]
seurat1[["umap.int"]]<-CreateDimReducObject(embeddings = umap1, key = "umap_int_", assay = "integrated")
seurat2[["umap.int"]]<-CreateDimReducObject(embeddings = umap2, key = "umap_int_", assay = "integrated")


# Save the original objects with the new coordinates
saveRDS(seurat1, file = seurat1.file)
saveRDS(seurat2, file = seurat2.file)
