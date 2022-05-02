##### Supp Figure 4 ###########
# Palantir UMAPs
  # Pseudotime
  # Differentiation potential
  # Branch probabilities

# Palantir violin plots 
##########################
library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(plyr)
source("/home/mainciburu/scRNA/colors.r")

# Young
seurat<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
plot_name<-"/home/mainciburu/scRNA/figures/supp_figure4/young_"
res_path<-"/home/mainciburu/scRNA/palantir/results/young/"

branch_prob<-read.csv(paste0(res_path, "branch_probs.csv"), row.names = 1)
pseudotime<-read.csv(paste0(res_path, "pseudotime.csv"), row.names = 1)
diff_pot<-read.csv(paste0(res_path, "diff_potential.csv"), row.names = 1)
colnames(diff_pot)<-"DP"
colnames(pseudotime)<-"pseudotime"


cols<-viridis::viridis(30)
reduction<-"umap.int"
# Pseudotime
pdf(file = paste0(plot_name, "pseudotime.pdf"), useDingbats = F,
    width = 6, height = 5)
seurat$plot<-pseudotime$pseudotime
FeaturePlot(seurat, reduction = reduction, features = "plot",
                  pt.size = 0.2) + ggtitle("Pseudotime") +
    scale_color_gradientn(colours = cols) + labs(colour = "Pseudotime")  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# Differentiation potential
pdf(file = paste0(plot_name, "diff_pot.pdf"), useDingbats = F,
    width = 6, height = 5)
seurat$plot<-diff_pot$DP
FeaturePlot(seurat, reduction = reduction, features = "plot",
            pt.size = 0.2) + ggtitle("Differentiation Potential") +
  scale_color_gradientn(colours = cols) + labs(colour = "Differentiation\nPotential")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# Branch probabilities
for(i in 1:ncol(branch_prob)){
  branch<-colnames(branch_prob)[i]
  pdf(file = paste0(plot_name, branch, "_branch_probability.pdf"), useDingbats = F,
      width = 6, height = 5)
  seurat$plot<-branch_prob[,i]
  print(FeaturePlot(seurat, reduction = reduction, features = "plot",
              pt.size = 0.2) + ggtitle(branch) +
    scale_color_gradientn(colours = cols) + labs(colour = paste0(branch, " \nBranch \nProbability"))  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2"))
  dev.off()
}

# Senior
seurat<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")
plot_name<-"/home/mainciburu/scRNA/figures/supp_figure4/senior_"
res_path<-"/home/mainciburu/scRNA/palantir/results/senior/"

branch_prob<-read.csv(paste0(res_path, "branch_probs.csv"), row.names = 1)
pseudotime<-read.csv(paste0(res_path, "pseudotime.csv"), row.names = 1)
diff_pot<-read.csv(paste0(res_path, "diff_potential.csv"), row.names = 1)
colnames(diff_pot)<-"DP"
colnames(pseudotime)<-"pseudotime"

cols<-viridis::viridis(30)
reduction<-"umap.int"
# Pseudotime
pdf(file = paste0(plot_name, "pseudotime.pdf"), useDingbats = F,
    width = 6, height = 5)
seurat$plot<-pseudotime$pseudotime
FeaturePlot(seurat, reduction = reduction, features = "plot",
                  pt.size = 0.2) + ggtitle("Pseudotime") +
    scale_color_gradientn(colours = cols) + labs(colour = "Pseudotime")  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# Differentiation potential
pdf(file = paste0(plot_name, "diff_pot.pdf"), useDingbats = F,
    width = 6, height = 5)
seurat$plot<-diff_pot$DP
FeaturePlot(seurat, reduction = reduction, features = "plot",
            pt.size = 0.2) + ggtitle("Differentiation Potential") +
  scale_color_gradientn(colours = cols) + labs(colour = "Differentiation\nPotential")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# Branch probabilities
for(i in 1:ncol(branch_prob)){
  branch<-colnames(branch_prob)[i]
  pdf(file = paste0(plot_name, branch, "_branch_probability.pdf"), useDingbats = F,
      width = 6, height = 5)
  seurat$plot<-branch_prob[,i]
  print(FeaturePlot(seurat, reduction = reduction, features = "plot",
              pt.size = 0.2) + ggtitle(branch) +
    scale_color_gradientn(colours = cols) + labs(colour = paste0(branch, " \nBranch \nProbability"))  + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
    theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2"))
  dev.off()
}

# Palantir violin plots - go to 04_trajectory_analysis/02_palantir/07_palantir_stats.R
