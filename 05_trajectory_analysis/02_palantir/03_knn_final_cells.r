library(Seurat)
# Final cells young
branch.prob.file<-"/home/mainciburu/scRNA/palantir/results/young/branch_probs.csv"
branch.prob<-read.csv(file = branch.prob.file, row.names = 1)
end.cells<-apply(branch.prob, 2, function(x){names(x)[which.max(x)]})

# umap coordinates for end cells in young
umap.file<-"/home/mainciburu/scRNA/palantir/data/young_umap.txt"
umap.coord<-read.table(umap.file, row.names = 1)

# seurat object condition 2
seurat.file<-"/home/mainciburu/scRNA/senior/seurat_senior_v4.rds"
seurat<-readRDS(seurat.file)

# query -> young cells in integrated umap
df<-umap.coord[end.cells,]
df<-as.data.frame(df)
df$CellType<-names(end.cells)

# data -> umap coordinates for cells in condition 2
condition2<-"senior"
dat<-seurat@reductions$umap.int@cell.embeddings
# results
res<-data.frame(query = rownames(df), query.ident=df$CellType,
                res = NA, res.ident=NA, distance=NA)
for(i in 1:nrow(df)){
  k<-FNN::get.knnx(data = dat, query = df[i,1:2], k = 5)
  s<-rownames(dat)[k$nn.index[1,]]
  names(s)<-as.character(seurat$prediction[s])
  si<-which(names(s)==df[i,3])[1]     # cell from same type
  si<-ifelse(is.na(si), yes = 1, no = si)   # if there is none -> chose 1st nn
  s<-s[si]
  res[i,3]<-s
  res[i,4]<-names(s)
  res[i,5]<-k$nn.dist[1,si]
}
print(res)
# Save cells as csv
res.file<-"/home/mainciburu/scRNA/palantir/data/final_cells_senior.csv"
write.csv(res$res, file = res.file, row.names = F, quote = F)
