######################## Erythroid ########################

library(pheatmap)
library(Seurat)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(ggalluvial)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(xlsx)
library(plyr)
library(grid)
source("/home/mainciburu/scRNA/colors.r")
source("/home/mainciburu/scRNA/pipeline_2.0/08_trajectory_analysis/plot_gene_trends.r")
source("/home/mainciburu/scRNA/pipeline_2.0/08_trajectory_analysis/compute_gene_trends.r")
source("/home/mainciburu/scRNA/ScoreSet.r")

pst.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/pseudotime.csv", row.names = 1)
pst.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/pseudotime.csv", row.names = 1)
pst.m<-read.csv("/home/mainciburu/scRNA/palantir/results/mds5/pseudotime.csv", row.names = 1)

bp.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/branch_probs.csv", row.names = 1)
bp.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/branch_probs.csv", row.names = 1)
bp.m<-read.csv("/home/mainciburu/scRNA/palantir/results/mds5/branch_probs.csv", row.names = 1)

res.y<-readRDS("scRNA/palantir/gene_trends/young/results_trends_young_ery.rds")
res.s<-readRDS("scRNA/palantir/gene_trends/senior/results_trends_senior_ery.rds")
res.m<-readRDS("scRNA/palantir/gene_trends/mds5/results_trends_mds5_ery.rds")

plot_path<-"/home/mainciburu/scRNA/palantir/pics/"
branch<-"Erythroid_late"

######## Select genes of interest
# markers of any cluster inside the trajectory
deg.young<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_young.rds")
deg.senior<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_senior.rds")
deg.mds<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_mds5.rds")

deg.young<-deg.young[deg.young$p_val_adj<=0.01,]
deg.senior<-deg.senior[deg.senior$p_val_adj<=0.01,]
deg.mds<-deg.mds[deg.mds$p_val_adj<=0.01,]

deg.ery<-deg.young$gene[deg.young$cluster%in%c("HSC", "MEP", "Erythroid_early", "Erythroid_late") & 
                           deg.young$avg_logFC>0.4]
deg.ery<-c(deg.ery, 
            deg.senior$gene[deg.senior$cluster%in%c("HSC", "MEP", "Erythroid_early", "Erythroid_late") & 
                           deg.senior$avg_logFC>0.4])
deg.ery<-c(deg.ery, 
            deg.mds$gene[deg.mds$cluster%in%c("HSC", "MEP", "Erythroid_early", "Erythroid_late") & 
                           deg.mds$avg_logFC>0.4])
deg.ery<-unique(deg.ery)
deg.ery<-intersect(deg.ery, rownames(res.y[[branch]]$trends))
deg.ery<-intersect(deg.ery, rownames(res.s[[branch]]$trends))
deg.ery<-intersect(deg.ery, rownames(res.m[[branch]]$trends))

######################################################
#################### Clusters ########################
cl.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/clusters_Erythroid_late.csv", 
               header = T, stringsAsFactors = F)
cl.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/clusters_Erythroid_late.csv", 
               header = T, stringsAsFactors = F)
cl.m<-read.csv("/home/mainciburu/scRNA/palantir/results/mds5/clusters_Erythroid_late.csv", 
               header = T, stringsAsFactors = F)
colnames(cl.y)<-colnames(cl.s)<-colnames(cl.m)<-c("V1", "V2")
cl.y<-cl.y[cl.y$V1%in%deg.ery,]
cl.s<-cl.s[cl.s$V1%in%deg.ery,]
cl.m<-cl.m[cl.m$V1%in%deg.ery,]

## order and merge similar clusters
cl.order<-c(8, 11, 9, 10, 2, 6, 7, 1, 0, 12, 4, 3, 5)
ii<-c()
for(cc in cl.order){ii<-c(ii, which(cl.y$V2==cc))}
cl.y<-cl.y[ii,]
cl.y$V2<-mapvalues(x = cl.y$V2, from = cl.order, to = 0:12)

x<-res.y[["Erythroid_late"]]$trends[cl.y$V1,]
col<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(30)
annCol<-c(RColorBrewer::brewer.pal(12, "Set3")[1:11], "#FF7F00", "#B15928")
df<-data.frame(cluster=factor(cl.y$V2), row.names = cl.y$V1)
annCol<-list(cluster=annCol)
names(annCol$cluster)<-as.character(0:12)
pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F, show_rownames = F,
         width = 6, height = 10, annotation_colors = annCol,
         filename = paste0(plot_path, "young_heatmap_trends_ery_deg.pdf"))

x<-res.s[["Erythroid_late"]]$trends[cl.y$V1,]
cl.order<-c(2, 8, 3, 9, 4, 5, 11, 10, 0, 1, 7, 6)
cl.s<-cl.s[match(cl.y$V1, cl.s$V1),]
cl.s$V2<-mapvalues(x = cl.s$V2, from = cl.order, to = 0:11)

df<-data.frame(cluster=factor(cl.s$V2), row.names = cl.s$V1)
pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F, show_rownames = F, 
         width = 6, height = 10, annotation_colors = annCol,
         filename = paste0(plot_path, "senior_heatmap_trends_ery_deg.pdf"))

x<-res.m[["Erythroid_late"]]$trends[cl.y$V1,]
cl.order<-c(5, 10, 8, 4, 0, 6, 11, 7, 2, 9, 3, 1)
cl.m<-cl.m[match(cl.y$V1, cl.m$V1),]
cl.m$V2<-mapvalues(x = cl.m$V2, from = cl.order, to = 0:11)

df<-data.frame(cluster=factor(cl.m$V2), row.names = cl.m$V1)
pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F, show_rownames = F, 
         width = 6, height = 10, annotation_colors = annCol,
         filename = paste0(plot_path, "mds_heatmap_trends_ery_deg.pdf"))



######################################################
### Correlater gene trends and branch probability trends
res.bp.y<-compute_gene_trends(branch_prob = bp.y, pseudotime = pst.y,
                             gene_exprs = bp.y, lineages = colnames(bp.y), ncores = 3,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/young/bp/")

res.bp.s<-compute_gene_trends(branch_prob = bp.s, pseudotime = pst.s,
                             gene_exprs = bp.s, lineages = colnames(bp.s), ncores = 3,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/senior/bp/")

res.bp.m<-compute_gene_trends(branch_prob = bp.m, pseudotime = pst.m,
                             gene_exprs = bp.m, lineages = colnames(bp.m), ncores = 3,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/mds5/bp/")
genes<-res.y[[branch]]$trends[deg.ery,]
bp.trend<-res.bp.y[[branch]]$trends[branch,]
cor.y<-apply(genes, 1, function(x){cor(x, bp.trend)})

genes<-res.s[[branch]]$trends[deg.ery,]
bp.trend<-res.bp.s[[branch]]$trends[branch,]
cor.s<-apply(genes, 1, function(x){cor(x, bp.trend)})

genes<-res.m[[branch]]$trends[deg.ery,]
bp.trend<-res.bp.m[[branch]]$trends[branch,]
cor.m<-apply(genes, 1, function(x){cor(x, bp.trend)})

df<-data.frame(gene=deg.ery, cor.young=cor.y, cor.s = cor.s, cor.m = cor.m,
                y_vs_s=abs(cor.y - cor.s), y_vs_m=abs(cor.y - cor.m), s_vs_m=abs(cor.s - cor.m))
df$cl.y<-cl.y$V2[match(deg.ery, cl.y$V1)]
df$cl.s<-cl.s$V2[match(deg.ery, cl.s$V1)]
df$cl.m<-cl.m$V2[match(deg.ery, cl.m$V1)]

# Select TFs
library(biomaRt)
mart<-useEnsembl("ensembl")
mart<-useDataset(dataset = "hsapiens_gene_ensembl", mart = mart)

c<-getBM(attributes = c("hgnc_symbol", "name_1006"), filters = "hgnc_symbol", values = deg.ery, mart = mart, verbose = T)
tf<-c$hgnc_symbol[c$name_1006%in%c("DNA-binding transcription factor activity",
                                   "DNA-binding transcription factor activity, RNA polymerase II-specific")]

df$TF<-0
df$TF[df$gene%in%tf]<-1

write.xlsx(df, "/home/mainciburu/scRNA/palantir/erythroid_trend_branchProb_correlation.xlsx", row.names=F)

### plot trends for every deg.ery in young, senior and MDS
pdf(paste0(plot_path, "gene_trends_ery/deg_ery.pdf"), width = 8, height = 6)
for(g in deg.ery){
  p<-plot_trends(r1 = res.y, r2=res.s, r3 = res.m, gg = g, bb="Erythroid_late", mode = 311, 
                 col1=col.condition[1], col2=col.condition[2], col3 = col.condition[6],
                 c1="Young", c2="Elderly", c3="MDS")
  print(p)
}
dev.off()


###### MDS heatmap adding labels to plotted genes
gg<-c("SLC25A6", "EPSTI1", "PHF6", "JUN", "LMO2", "LYL1", "YBX1", "KLF1", "GATA1", "AHSP", "CA1", "HBB")
i<-match(cl.y$V1, gg)
i<-i[!is.na(i)]
gg<-gg[i]

x<-res.m[["Erythroid_late"]]$trends[cl.y$V1,]
annCol<-c(RColorBrewer::brewer.pal(12, "Set3")[1:11], "#FF7F00", "#B15928")
annCol<-list(cluster=annCol, plot=c("No"="grey", "Yes"="red"))
names(annCol$cluster)<-names(annCol$cluster)<-as.character(0:12)

df<-data.frame(cluster=factor(cl.m$V2), plot="No", row.names = cl.m$V1)
i<-match(gg, cl.m$V1)
df$plot<-as.character(df$plot)
df$plot[i]<-"Yes"

ph<-pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F, show_rownames = T, 
         annotation_colors = annCol)


source("/home/mainciburu/scRNA/palantir/add.flag.R")
ph<-add.flag(ph,
         kept.labels = gg,
         repel.degree = 0.5)

pdf(paste0(plot_path, "mds5_heatmap_trends_ery_deg_annot.pdf"),
    height = 8, width = 6)
grid.newpage()
grid.draw(ph)
dev.off()

