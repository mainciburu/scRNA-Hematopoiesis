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
source("/home/mainciburu/scRNA/pipeline_2.0/functions.r")
source("/home/mainciburu/scRNA/pipeline_2.0/04_trajectory_analysis/02_palantir/compute_gene_trends.r")

### Input
res_path<-"/home/mainciburu/scRNA/palantir/results/"
plot_path<-"/home/mainciburu/scRNA/palantir/pics/"
res.y<-readRDS("/home/mainciburu/scRNA/palantir/gene_trends/young/results_trends_young_Mono_pDC.rds")
res.s<-readRDS("/home/mainciburu/scRNA/palantir/gene_trends/senior/results_trends_senior_Mono_pDC.rds")

####### Select genes of interest
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


### Gene trend heatmaps
cl.y<-read.csv(paste0(res_path, "young/clusters_Monocytes.csv"), 
             header = T, stringsAsFactors = F)
colnames(cl.y)<-c("V1", "V2")
cl.y<-cl.y[cl.y$V1%in%deg.mono,]
g<-cl.y$V1

# Join similar clusters
cl.y$V3<-cl.y$V2
cl.y$V3[cl.y$V2%in%c(5, 7)]<-1
cl.y$V3[cl.y$V2==9]<-2
cl.y$V3[cl.y$V2%in%c(0, 3)]<-3
cl.y$V3[cl.y$V2%in%c(6, 4)]<-4
cl.y$V3[cl.y$V2%in%c(8, 10)]<-5
cl.y$V3[cl.y$V2==2]<-6
cl.y$V3[cl.y$V2==1]<-7
cl.y<-cl.y[order(cl.y$V3),]
g<-cl.y$V1

x<-res.y[["Monocytes"]]$trends[g,]
x<-x[g,]
annCol<-RColorBrewer::brewer.pal(length(unique(cl.y$V3)), "Set3")
df<-data.frame(cluster=factor(cl.y$V3), row.names = cl.y$V1)
annCol<-list(cluster=annCol)
names(annCol$cluster)<-as.character(1:length(unique(cl.y$V3)))

# set scale
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}
sc<-c(scale_rows(res.y[["Monocytes"]]$trends[cl.y$V1,]),
      scale_rows(res.s[["Monocytes"]]$trends[cl.y$V1,]))
Breaks<-seq(min(sc), max(sc), length = 100)
col<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df,
         cluster_cols = F, show_colnames = F, breaks = Breaks,
         show_rownames = F, width = 6.5, height = 10, annotation_colors = annCol,
         filename = paste0(plot_path, "young_heatmap_trends_monocytes_deg_merged.png"))
g.y<-g

# Senior
x<-res.s[["Monocytes"]]$trends[cl.y$V1,]
pheatmap(x, color = col, scale = "row", cluster_rows = F,
         cluster_cols = F, show_colnames = F, breaks = Breaks,
         show_rownames = F, width = 6, height = 10, annotation_colors = annCol,
         filename = paste0(plot_path, "senior_heatmap_trends_monocytes_deg_merged_young_order.png"))


###### Enrichment Young per merged cluster ######
clust<-cl.y[,c(1,3)]
colnames(clust)<-c("gene", "cluster")

# get gene sets
m_df = msigdbr(species = "Homo sapiens")
m_df = m_df[m_df$gs_cat == "H" | m_df$gs_subcat %in% c("BP", "CP", "CP:REACTOME", "CP:BIOCARTA", "CP:PID", "CP:KEGG"),]
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# background genes
seurat<-readRDS("scRNA/young/seurat_young_v3.rds")
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
write.xlsx(x = df, file = "/home/mainciburu/scRNA/palantir/Monocytes_enrichment_degMono_merged_young.xlsx")


#### Correlate gene trend and branch probability 
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

genes<-res.y[["Monocytes"]]$trends[g.y,]
bp.trend<-res.bp.y[["Monocytes"]]$trends["Monocytes",]
cor.y<-apply(genes, 1, function(x){cor(x, bp.trend)})

genes<-res.s[["Monocytes"]]$trends[g.y,]
bp.trend<-res.bp.s[["Monocytes"]]$trends["Monocytes",]
cor.s<-apply(genes, 1, function(x){cor(x, bp.trend)})

df<-data.frame(gene=g.y, cor.young=cor.y, cor.s = cor.s, 
                difference=abs(cor.y - cor.s), cl.young=cl.y[g.y,3])
df$TF<-0
df$TF[df$gene%in%tf]<-1

write.xlsx(df, "/home/mainciburu/scRNA/palantir/monocytes_trend_branchProb_correlation.xlsx", row.names=F)

# plot trends for every deg.mono in young and senior
pdf(paste0(plot_path, "gene_trends_mono/deg_mono.pdf"), width = 6, height = 10)
for(gi in g.y){
    df.y<-data.frame(trend=res.y$Monocytes$trends[gi,],
                 std=res.y$Monocytes$std[gi,],
                 condition="young")
    df.y$ymin<-df.y$trend - df.y$std
    df.y$ymax<-df.y$trend + df.y$std
    df.y$time<-as.numeric(rownames(df.y))

    df.s<-data.frame(trend=res.s$Monocytes$trends[gi,],
                 std=res.s$Monocytes$std[gi,],
                 condition="senior")
    df.s$ymin<-df.s$trend - df.s$std
    df.s$ymax<-df.s$trend + df.s$std
    df.s$time<-as.numeric(rownames(df.s))

    df<-rbind(df.y, df.s)

    print(plot_trend_2conditions(df, gi))
}


# plot trends for selected genes and stats table

gg<-c("CST7", "PRTN3", "MPO", "CD74", "CALR", "GNAS", "MYCT1", "MLLT3", "ALDH1A1", "FOS", "JUNB", "TSC22D3", "DUSP1", "DDIT4", "ANXA1", "LGALS3", "JAML")
stats.res<-matrix(NA, nrow = length(g.y), ncol = 3)
rownames(stats.res)<-g.y
colnames(stats.res)<-c("Wilcox_PValue", "Wilcox_Stat", "ToTalVariation")

for(gi in g.y){
    df.y<-data.frame(trend=res.y$Monocytes$trends[gi,],
                 std=res.y$Monocytes$std[gi,],
                 condition="young")
    df.y$ymin<-df.y$trend - df.y$std
    df.y$ymax<-df.y$trend + df.y$std
    df.y$time<-as.numeric(rownames(df.y))

    df.s<-data.frame(trend=res.s$Monocytes$trends[gi,],
                 std=res.s$Monocytes$std[gi,],
                 condition="senior")
    df.s$ymin<-df.s$trend - df.s$std
    df.s$ymax<-df.s$trend + df.s$std
    df.s$time<-as.numeric(rownames(df.s))
    df<-rbind(df.y, df.s)
    test_results <- wilcox.test(df$trend ~ df$condition)
    tv_value <- sum(abs(df$trend[df$condition=="young"] - df$trend[df$condition=="senior"]))
    stats.res[gi,1]<-test_results$p.val
    stats.res[gi,2]<-test_results$statistic
    stats.res[gi,3]<-tv_value
    #print(plot_trend_2conditions(df, gi))
}
stats.res<-data.frame(stats.res)
stats.res$Wilcox_PValue_adjusted<-p.adjust(stats.res$Wilcox_PValue, method = "BH")
write.csv(stats.res, file = "/home/mainciburu/scRNA/palantir/monocytes_gene_trends_stats.csv", quote = F)

### senior heatmap labeling plotted genes
x<-res.s[["Monocytes"]]$trends[cl.y$V1,]
df<-data.frame( plot="No", row.names = rownames(x))
annCol<-list(cluster=annCol, plot=c("No"="grey", "Yes"="red"))
i<-match(gg, rownames(x))
df$plot<-as.character(df$plot)
df$plot[i]<-"Yes"

ph<-pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df, breaks = Breaks,
         cluster_cols = F,  show_colnames = F, annotation_colors = annCol)

source("/home/mainciburu/scRNA/palantir/add.flag.R")

i<-match(rownames(x), gg)
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
