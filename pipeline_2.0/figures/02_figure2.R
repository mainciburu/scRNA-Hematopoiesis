##### Figure 2 ###########
# STREAM plots
# Pseudotime vs branch probability dotplot
# Gene trends heatmap
# Gene trends lineplot

##########################
library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(plyr)
library(pheatmap)

source("/home/mainciburu/scRNA/colors.r")
young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")

### STREAM plots - go to 04_trajectory_analysis/01_stream/Stream.py

#### Pseudotime vs branch probability dotplot
# Young
res_path<-"/home/mainciburu/scRNA/palantir/results/young/"
plot_path<-"/home/mainciburu/scRNA/figures/figure2/"
branch_prob<-read.csv(paste0(res_path, "branch_probs.csv"), row.names = 1)
pseudotime<-read.csv(paste0(res_path, "pseudotime.csv"), row.names = 1)
colnames(pseudotime)<-"pseudotime"
df<-data.frame(pst=pseudotime$pseudotime,
               group = young$CellType2,
               patient=young$Patient)
df<-cbind(df, Branch = branch_prob[,"Monocytes"])
pp<-ggplot(df, aes(pst, Branch, colour = group)) + geom_point(size = 0.7) + theme_bw() +
    scale_color_manual(values = col.young.v3) + labs(colour = "Cell Type")  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "Pseudotime", y = "Branch Probability") +
    ggtitle(branch)
pdf(paste0(plot_path, "young_pst_vs_monocytes_branch_prob.pdf"), useDingbats = F, width = 12)
print(pp)
dev.off()

# Senior
res_path<-"/home/mainciburu/scRNA/palantir/results/senior/"
branch_prob<-read.csv(paste0(res_path, "branch_probs.csv"), row.names = 1)
pseudotime<-read.csv(paste0(res_path, "pseudotime.csv"), row.names = 1)
colnames(pseudotime)<-"pseudotime"
df<-data.frame(pst=pseudotime$pseudotime,
               group = senior$prediction,
               patient=senior$Patient)
df<-cbind(df, Branch = branch_prob[,"Monocytes"])
pp<-ggplot(df, aes(pst, Branch, colour = group)) + geom_point(size = 0.7) + theme_bw() +
    scale_color_manual(values = col.young.v3) + labs(colour = "Cell Type")  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "Pseudotime", y = "Branch Probability") +
    ggtitle(branch)
pdf(paste0(plot_path, "senior_pst_vs_monocytes_branch_prob.pdf"), useDingbats = F, width = 12)
print(pp)
dev.off()


### Gene trends heatmap: go to 04_trajectory_analysis/02_palantir/07_monocyte_analysis.r

### Gene trends lineplot: go to 04_trajectory_analysis/02_palantir/07_monocyte_analysis.r

