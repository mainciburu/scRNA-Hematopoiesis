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
library(biomaRt)
source("/home/mainciburu/scRNA/colors.r")
source("/home/mainciburu/scRNA/pipeline_2.0/functions.r")
source("/home/mainciburu/scRNA/pipeline_2.0/04_trajectory_analysis/02_palantir/compute_gene_trends.r")


######## Select genes of interest
# markers of any cluster inside the trajectory
deg.young<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_young.rds")
deg.senior<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_senior.rds")
deg.m1<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_mds1.rds")
deg.m2<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_mds2.rds")
deg.m3<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_mds3.rds")
deg.m4<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_mds4.rds")

res.y<-readRDS("scRNA/palantir/gene_trends/young/results_trends_young_ery.rds")
res.s<-readRDS("scRNA/palantir/gene_trends/senior/results_trends_senior_ery.rds")
res.m1<-readRDS("scRNA/palantir/gene_trends/mds1/results_trends_mds1_ery.rds")
res.m2<-readRDS("scRNA/palantir/gene_trends/mds2/results_trends_mds2_ery.rds")
res.m3<-readRDS("scRNA/palantir/gene_trends/mds3/results_trends_mds3_ery.rds")
res.m4<-readRDS("scRNA/palantir/gene_trends/mds4/results_trends_mds4_ery.rds")
res<-list(res.y,res.s,res.m1,res.m2,res.m3,res.m4)
names(res)<-c("Young", "Elderly", "MDS1", "MDS2", "MDS3", "MDS4")

deg.ery<-c()
for(deg in list(deg.young, deg.senior, deg.m1, deg.m2, deg.m3, deg.m4)){
   x<-deg[deg$p_val_adj<=0.01,]
   deg.ery<-c(deg.ery, x$gene[x$cluster%in%c("HSC", "MEP", "Erythroid_early", "Erythroid_late") & 
                           x$avg_logFC>0.4])
}

deg.ery<-unique(deg.ery)
deg.ery<-deg.ery[deg.ery%in%rownames(res.y[["Erythroid_late"]]$trends)]
deg.ery<-deg.ery[deg.ery%in%rownames(res.s[["Erythroid_late"]]$trends)]
deg.ery<-deg.ery[deg.ery%in%rownames(res.m1[["Erythroid_late"]]$trends)]
deg.ery<-deg.ery[deg.ery%in%rownames(res.m2[["Erythroid_late"]]$trends)]
deg.ery<-deg.ery[deg.ery%in%rownames(res.m3[["Erythroid_late"]]$trends)]
deg.ery<-deg.ery[deg.ery%in%rownames(res.m4[["Erythroid_late"]]$trends)]


### Analysis per patient
patients<-c("mds1", "mds2", "mds3", "mds4")

pst.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/pseudotime.csv", row.names = 1)
pst.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/pseudotime.csv", row.names = 1)

bp.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/branch_probs.csv", row.names = 1)
bp.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/branch_probs.csv", row.names = 1)

plot_path<-"/home/mainciburu/scRNA/MDS_paper/plots/trajectories/"
branch<-"Erythroid_late"

##################################################################################################
#################### Heatmaps ########################
cl.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/clusters_Erythroid_late.csv", 
               header = T, stringsAsFactors = F)

colnames(cl.y)<-c("V1", "V2")
## order and merge similar clusters
cl.order<-c(8, 11, 9, 10, 2, 6, 7, 1, 0, 12, 4, 3, 5)
ii<-c()
for(cc in cl.order){ii<-c(ii, which(cl.y$V2==cc))}
cl.y<-cl.y[ii,]
cl.y$V2<-mapvalues(x = cl.y$V2, from = cl.order, to = 0:12)

cl.y<-cl.y[cl.y$V1%in%deg.ery,]

# set scale
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}
sc<-c(scale_rows(res.y[["Erythroid_late"]]$trends[cl.y$V1,]),
      scale_rows(res.s[["Erythroid_late"]]$trends[cl.y$V1,]),
      scale_rows(res.m1[["Erythroid_late"]]$trends[cl.y$V1,]),
      scale_rows(res.m2[["Erythroid_late"]]$trends[cl.y$V1,]),
      scale_rows(res.m3[["Erythroid_late"]]$trends[cl.y$V1,]),
      scale_rows(res.m4[["Erythroid_late"]]$trends[cl.y$V1,]))
Breaks<-seq(min(sc), max(sc), length = 100)
col<-colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

for(ix in 1:length(res)){
   if(ix==1){    # Young
      x<-res.y[["Erythroid_late"]]$trends[cl.y$V1,]
      annCol<-c(RColorBrewer::brewer.pal(12, "Set3")[1:11], "#FF7F00", "#B15928")
      df<-data.frame(cluster=factor(cl.y$V2), row.names = cl.y$V1)
      annCol<-list(cluster=annCol)
      names(annCol$cluster)<-as.character(0:12)
      pheatmap(x, color = col, scale = "row", cluster_rows = F, annotation_row = df,
               cluster_cols = F, show_colnames = F, show_rownames = F,
               width = 6, height = 10, annotation_colors = annCol, breaks = Breaks,
               filename = paste0(plot_path, names(res)[ix], "_heatmap_trends_ery_deg.png"))

   }else{
   x<-res[[ix]][["Erythroid_late"]]$trends[cl.y$V1,]
   pheatmap(x, color = col, scale = "row", cluster_rows = F, 
            cluster_cols = F, show_colnames = F, show_rownames = F, 
            width = 6, height = 10, breaks = Breaks,
            filename = paste0(plot_path, names(res)[ix], "_heatmap_trends_ery_deg.png"))
   }
}


######################################################
### Correlater gene trends and branch probability trends
res.bp.y<-compute_gene_trends(branch_prob = bp.y, pseudotime = pst.y,
                             gene_exprs = bp.y, lineages = "Erythroid_late", ncores = 3,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/young/bp/")

res.bp.s<-compute_gene_trends(branch_prob = bp.s, pseudotime = pst.s,
                             gene_exprs = bp.s, lineages = "Erythroid_late", ncores = 3,
                             res_path = "/home/mainciburu/scRNA/palantir/gene_trends/senior/bp/")


for(patient in patients){
   res.m<-readRDS(paste0("scRNA/palantir/gene_trends/", patient, "/results_trends_", patient, "_ery.rds"))
   pst.m<-read.csv(paste0("/home/mainciburu/scRNA/palantir/results/", patient, "/pseudotime.csv"), row.names = 1)
   bp.m<-read.csv(paste0("/home/mainciburu/scRNA/palantir/results/", patient, "/branch_probs.csv"), row.names = 1)
   res.bp.m<-compute_gene_trends(branch_prob = bp.m, pseudotime = pst.m,
                                gene_exprs = bp.m, lineages = "Erythroid_late", ncores = 3,
                                res_path = paste0("/home/mainciburu/scRNA/palantir/gene_trends/", patient, "/bp/"))
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

   # Select TFs
   mart<-useEnsembl("ensembl")
   mart<-useDataset(dataset = "hsapiens_gene_ensembl", mart = mart)

   c<-getBM(attributes = c("hgnc_symbol", "name_1006"), filters = "hgnc_symbol", values = deg.ery, mart = mart, verbose = T)
   tf<-c$hgnc_symbol[c$name_1006%in%c("DNA-binding transcription factor activity",
                                      "DNA-binding transcription factor activity, RNA polymerase II-specific")]
   df$TF<-0
   df$TF[df$gene%in%tf]<-1

   write.xlsx(df, paste0("/home/mainciburu/scRNA/palantir/erythroid_trend_branchProb_correlation_", patient, ".xlsx"), row.names=F)
}


### plot trends for every deg.ery in young, senior and MDS
pdf(paste0(plot_path, "trends_deg_ery.pdf"), width = 8, height = 6)
for(gi in deg.ery){
   df<-data.frame()
   for(i in 1:length(res)){
      dfi<-data.frame(trend=res[[i]]$Erythroid_late$trends[gi,],
                    std=res[[i]]$Erythroid_late$std[gi,],
                    condition=names(res)[i])
       dfi$ymin<-dfi$trend - dfi$std
       dfi$ymax<-dfi$trend + dfi$std
       dfi$time<-as.numeric(rownames(dfi))
       df<-rbind(df, dfi)
   }    
    print(plot_trend_multiple_conditions(df, gi))
}
dev.off()

### Statistic test
conditions<-c("Young", "Elderly", "MDS1", "MDS2", "MDS3", "MDS4")
condition.comb<-combn(conditions, 2)
wilcox.res<-matrix(NA, ncol = ncol(condition.comb), nrow = length(deg.ery))
rownames(wilcox.res)<-deg.ery
colnames(wilcox.res)<-paste0(condition.comb[1,],"_vs_",condition.comb[2,])
wilcox.pval<-wilcox.res
total.var<-wilcox.pval
for(gi in deg.ery){
   for(jj in 1:ncol(condition.comb)){
      c1<-condition.comb[1,jj]
      c2<-condition.comb[2,jj]
      df1<-data.frame(trend=res[[c1]]$Erythroid_late$trends[gi,],
                    condition=c1)
      df2<-data.frame(trend=res[[c2]]$Erythroid_late$trends[gi,],
                    condition=c2)
      df<-rbind(df1, df2)
      test_results <- wilcox.test(df$trend ~ df$condition)
      tv_value <- sum(abs(df$trend[df$condition==c1] - df$trend[df$condition==c2]))
      wilcox.pval[gi,jj]<-test_results$p.val
      wilcox.res[gi,jj]<-test_results$statistic
      total.var[gi,jj]<-tv_value
   }
}
wilcox.pval.adjust<-matrix(p.adjust(as.vector(wilcox.pval), method='BH'),ncol=ncol(wilcox.pval))
stats.res<-data.frame(Genes=rownames(stats.res))
for(ii in 1:ncol(wilcox.res)){
   dfi<-cbind(wilcox.res[,ii], wilcox.pval[,ii], wilcox.pval.adjust[,ii], total.var[,ii])
   colnames(dfi)<-paste0(c("Wilcox_Stat_", "Wilcox_PValue_", "Wilcox_PValue_adjusted_", "TotalVariation_"), colnames(wilcox.res)[ii])
   stats.res<-cbind(stats.res, dfi)
}
write.csv(stats.res, file = "scRNA/palantir/erythroid_gene_trends_stats.csv")

###### MDS 4 heatmap adding labels to plotted genes
gg<-c("YBX1", "HMGA1", "HK1", "JUN", "PHF6", "PDCD4", "TRIB2", "PVT1", "HBD", "GATA1", "AHSP")
i<-match(cl.y$V1, gg)
i<-i[!is.na(i)]
gg<-gg[i]

x<-res.m4[["Erythroid_late"]]$trends[cl.y$V1,]
annCol<-list(plot=c("No"="grey", "Yes"="red"))

df<-data.frame(plot=rep("No", length(cl.y$V1)), row.names = cl.y$V1)
i<-match(gg, cl.y$V1)
df$plot<-as.character(df$plot)
df$plot[i]<-"Yes"

ph<-pheatmap(x, color = col, scale = "row", cluster_rows = F, 
         annotation_row = df, breaks = Breaks,
         cluster_cols = F, show_colnames = F, show_rownames = T, 
         annotation_colors = annCol)


source("/home/mainciburu/scRNA/palantir/add.flag.R")
ph<-add.flag(ph,
         kept.labels = gg,
         repel.degree = 0.5)

pdf(paste0(plot_path, "mds4_heatmap_trends_ery_deg_annot.pdf"),
    height = 8, width = 6)
grid.newpage()
grid.draw(ph)
dev.off()

