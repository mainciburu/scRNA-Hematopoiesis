########### Monocytes ################

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
library(grid)
source("/home/mainciburu/scRNA/colors.r")
source("/home/mainciburu/scRNA/pipeline_2.0/07_trajectory_analysis/plot_gene_trends.r")
source("/home/mainciburu/scRNA/ScoreSet.r")
source("/home/mainciburu/scRNA/pipeline_2.0/07_trajectory_analysis/compute_gene_trends.r")

###################### Young ###################################
### Input
branch<-"Monocytes"
condition<-"young"
res_path<-"/home/mainciburu/scRNA/palantir/results/young/"
plot_name<-"/home/mainciburu/scRNA/palantir/pics/young_"

cl<-read.csv(paste0(res_path, "clusters_", branch,".csv"), 
             header = T, stringsAsFactors = F)
colnames(cl)<-c("V1", "V2")
res<-readRDS(paste0("/home/mainciburu/scRNA/palantir/gene_trends/", 
                    condition, "/results_trends_", condition, "_Mono_pDC.rds"))
res<-res[[branch]]$trends

# Select genes of interest
# markers of any cluster inside the trajectory
deg.young<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_young.rds")
deg.senior<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_senior.rds")
deg.young<-deg.young[deg.young$p_val_adj<=0.01,]
deg.senior<-deg.senior[deg.senior$p_val_adj<=0.01,]

deg.mono<-deg.young$gene[deg.young$cluster%in%c("HSC", "LMPP", "GMP", "GMP_Granulocytes", "Monocytes") & 
                           deg.young$avg_logFC>0.4]
deg.mono<-c(deg.mono, 
            deg.senior$gene[deg.senior$cluster%in%c("HSC", "LMPP", "GMP", "GMP_Granulocytes", "Monocytes") & 
                           deg.senior$avg_logFC>0.4])
deg.mono<-unique(deg.mono)


cl<-cl[cl$V1%in%deg.mono,]
g<-cl$V1
res<-res[g,]

# Join similar clusters
cl$V3<-cl$V2
cl$V3[cl$V2%in%c(5, 7)]<-1
cl$V3[cl$V2==9]<-2
cl$V3[cl$V2%in%c(0, 3)]<-3
cl$V3[cl$V2%in%c(6, 4)]<-4
cl$V3[cl$V2%in%c(8, 10)]<-5
cl$V3[cl$V2==2]<-6
cl$V3[cl$V2==1]<-7
cl<-cl[order(cl$V3),]
g<-cl$V1
res<-res[g,]
annCol<-RColorBrewer::brewer.pal(length(unique(cl$V3)), "Set3")
df<-data.frame(cluster=factor(cl$V3), row.names = cl$V1)
annCol<-list(cluster=annCol)
names(annCol$cluster)<-as.character(1:length(unique(cl$V3)))
col<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(30)
pheatmap(res, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F,
         show_rownames = F, width = 5, height = 8, annotation_colors = annCol,
         filename = paste0(plot_name, "heatmap_trends_monocytes_deg_merged1.pdf"))
cl.y<-cl      # save for later
res.y<-res
g.y<-g


# Enrichment per merged cluster
clust<-cl[,c(1,3)]
colnames(clust)<-c("gene", "cluster")

# get gene sets
m_df = msigdbr(species = "Homo sapiens")
m_df = m_df[m_df$gs_cat == "H" | m_df$gs_subcat %in% c("BP", "CP", "CP:REACTOME", "CP:BIOCARTA", "CP:PID", "CP:KEGG"),]
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# background genes
bg<-rownames(seurat)

enrichres<-sapply(1:max(clust$cluster), function(x){
  clusterProfiler::enricher(gene = clust$gene[clust$cluster==x], universe = bg, TERM2GENE = m_t2g)
})
names(enrichres)<-1:(length(enrichres))
df<-data.frame()
for(i in 1:length(enrichres)){
  enrichres[[i]]@result$cluster<-names(enrichres)[i]
  df.i<-enrichres[[i]]@result
  df.i<-df.i[df.i$p.adjust<0.05,]
  df<-rbind(df, df.i)
}
write.xlsx(x = df, file = paste0("/home/mainciburu/scRNA/palantir/", branch, "_enrichment_degMono_merged_", condition, ".xlsx"))


###################### Senior ###################################
### Input
branch<-"Monocytes"
condition<-"senior"
res_path<-"/home/mainciburu/scRNA/palantir/results/senior/"
plot_name<-"/home/mainciburu/scRNA/palantir/pics/senior_"

cl<-read.csv(paste0(res_path, "clusters_", branch,".csv"), 
             header = T, stringsAsFactors = F)
colnames(cl)<-c("V1", "V2")
res<-readRDS(paste0("/home/mainciburu/scRNA/palantir/gene_trends/", 
                    condition, "/results_trends_", condition, "_Mono_pDC.rds"))

# Heatmap using Palantir gene clusters
res<-res[[branch]]$trend
cl<-cl[cl$V1%in%deg.mono,]

# Join similar clusters
cl$V3<-cl$V2
cl$V3[cl$V2==4]<-1
cl$V3[cl$V2==7]<-2
cl$V3[cl$V2%in%c(6, 2)]<-3
cl$V3[cl$V2==8]<-4
cl$V3[cl$V2%in%c(10, 0, 5, 3)]<-5
cl$V3[cl$V2==1]<-6
cl$V3[cl$V2==9]<-7

# heatmap with young clustering order
cl<-cl[match(g.y, cl$V1),]
res<-res[cl$V1,]
annCol<-RColorBrewer::brewer.pal(length(unique(cl$V3)), "Set3")
df<-data.frame(cluster=factor(cl$V3), row.names = cl$V1)
annCol<-list(cluster=annCol)
names(annCol$cluster)<-as.character(1:length(unique(cl$V3)))
pheatmap(res, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F,
         show_rownames = F, width = 8, height = 6, annotation_colors = annCol,
         filename = paste0(plot_name, "heatmap_trends_monocytes_deg_merged_young_order.pdf"))
cl.s<-cl    # save for later
res.s<-res


####################### Compare #############################
plot_path<-"/home/mainciburu/scRNA/palantir/pics/"

### Correlate gene trend and branch probability 
pst.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/pseudotime.csv", row.names = 1)
pst.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/pseudotime.csv", row.names = 1)

bp.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/branch_probs.csv", row.names = 1)
bp.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/branch_probs.csv", row.names = 1)

res.bp.y<-compute_gene_trends(branch_prob = bp.y, pseudotime = pst.y,
                             gene_exprs = bp.y, lineages = colnames(bp.y), ncores = 2,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/young/bp/")

res.bp.s<-compute_gene_trends(branch_prob = bp.s, pseudotime = pst.s,
                             gene_exprs = bp.s, lineages = colnames(bp.s), ncores = 2,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/senior/bp/")
rownames(cl.y)<-cl.y$V1
rownames(cl.s)<-cl.s$V1

branch<-"Monocytes"
genes<-res1[[branch]]$trends[g.y,]
bp.trend<-res.bp.y[[branch]]$trends[branch,]
cor.y<-apply(genes, 1, function(x){cor(x, bp.trend)})

genes<-res2[[branch]]$trends[g.y,]
bp.trend<-res.bp.s[[branch]]$trends[branch,]
cor.s<-apply(genes, 1, function(x){cor(x, bp.trend)})

df<-data.frame(gene=g.y, cor.young=cor.y, cor.s = cor.s, 
                difference=abs(cor.y - cor.s), cl.young=cl.y[g.y,3], cl.senior=cl.s[g.y,3])
df$TF<-0
df$TF[df$gene%in%tf]<-1

write.xlsx(df, "/home/mainciburu/scRNA/palantir/monocytes_trend_branchProb_correlation.xlsx", row.names=F)

# plot trends for every deg.mono in young and senior
res1<-readRDS(paste0("/home/mainciburu/scRNA/palantir/gene_trends/young/results_trends_young_Mono_pDC.rds"))
res2<-readRDS(paste0("/home/mainciburu/scRNA/palantir/gene_trends/senior/results_trends_senior_Mono_pDC.rds"))

pdf(paste0(plot_path, "gene_trends_mono/deg_mono.pdf"), width = 8, height = 6)
for(gi in g.y){
  p<-plot_trends(r1 = res1, r2=res2, gg = gi, bb="Monocytes", mode = 211, 
                 col1=col.condition[1], col2=col.condition[2],
                 c1="Young", c2="Elderly")
  print(p)
}
dev.off()

# plot trends for selected genes
gg<-c("CST7", "PRTN3", "MPO", "CD74", "CALR", "GNAS", "MYCT1", "MLLT3", "ALDH1A1", "FOS", "JUNB", "TSC22D3", "DUSP1", "DDIT4", "ANXA1", "LGALS3", "JAML")
for(gi in gg){
  pdf(paste0(plot_path, "gene_trends_mono/", gi, ".pdf"), width = 8, height = 6)
  p<-plot_trends(r1 = res1, r2=res2, gg = gi, bb="Monocytes", mode = 211, 
                 col1=col.condition[1], col2=col.condition[2],
                 c1="Young", c2="Elderly")
  print(p)
  dev.off()
}

### senior heatmap labeling plotted genes
cl.s<-cl.s[match(g.y, cl.s$V1),]
res.s<-res.s[cl.s$V1,]
annCol<-RColorBrewer::brewer.pal(length(unique(cl.s$V3)), "Set3")
df<-data.frame(cluster=factor(cl.s$V3), plot="No", row.names = cl.s$V1)
annCol<-list(cluster=annCol, plot=c("No"="grey", "Yes"="red"))
names(annCol$cluster)<-as.character(1:length(unique(cl.s$V3)))
i<-match(gg, cl.s$V1)
df$plot<-as.character(df$plot)
df$plot[i]<-"Yes"

ph<-pheatmap(res.s, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F,  show_colnames = F, annotation_colors = annCol)

source("/home/mainciburu/scRNA/palantir/add.flag.R")

i<-match(cl.s$V1, gg)
i<-i[!is.na(i)]
gg<-gg[i]

ph<-add.flag(ph,
         kept.labels = gg,
         repel.degree = 0.5)

pdf(paste0(plot_path, "senior_heatmap_trends_monocytes_deg_merged_young_order_annot.pdf"),
    height = 8, width = 6)
grid.newpage()
grid.draw(ph)
dev.off()
