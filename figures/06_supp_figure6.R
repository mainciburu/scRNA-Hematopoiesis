##### Supp Figure 6 ###########
# GSEA Young vs MDS
# UMAP erythroid branch probability in MDS
# Heatmaps gene trends
# Networks
##########################
library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(plyr)
source("/home/mainciburu/scRNA/colors.r")

# GSEA - go to 03_differential_expression/03_GSEA_young_elderly_mds.r

# UMAP erythroid branch probability in MDS
patients<-c("mds1", "mds2", "mds3", "mds4")
for(patient in patients){
	seurat<-readRDS(paste0("/home/mainciburu/scRNA/MDS_paper/seurat_", patient, ".rds"))
	plot_name<-paste0("/home/mainciburu/scRNA/figures/supp_figure4/", patient, "_")
	res_path<-paste0("/home/mainciburu/scRNA/palantir/results/", patient, "/")
	branch_prob<-read.csv(paste0(res_path, "branch_probs.csv"), row.names = 1)

	branch<-"Erythroid_late"
	pdf(file = paste0(plot_name, branch, "_branch_probability.pdf"), useDingbats = F,
	    width = 6, height = 5)
	seurat$plot<-branch_prob[,branch]
	print(FeaturePlot(seurat, reduction = reduction, features = "plot",
	            pt.size = 0.2) + ggtitle(branch) +
	  scale_color_gradientn(colours = cols) + labs(colour = paste0(branch, " \nBranch \nProbability"))  + 
	  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
	  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
	  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2"))
	dev.off()
}


# Heatmaps gene trends - go to 04_trajectory_analysis/02_palantir/09_erythroid_analysis.r

# Networks - go to 05_GRN/05_CytoscapeVisualization.R



