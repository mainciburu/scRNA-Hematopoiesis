######################
# GSEA per cluster between young, elderly and MDS
# MsigDb Hallmark pathways
# Plot results
##########################

library(fgsea)
library(ggplot2)
library(msigdbr)
library(foreach)
library(doParallel)
source("/home/mainciburu/scRNA/colors.r")

### Hallmark genesets
m_df = msigdbr(species = "Homo sapiens")
m_df = m_df[m_df$gs_cat == "H",]
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### GSEA per MDS patient
patients<-c("MDS1", "MDS2", "MDS3", "MDS4")

for(patient in patients){

  ### Input
  deg.file<-paste0("/home/mainciburu/scRNA/btwn_condition/deg_MAST_latent_var_young_", patient, ".Rdata")
  load(deg.file)
  deg1<-deg.condition
  deg.file<-paste0("/home/mainciburu/scRNA/btwn_condition/deg_MAST_latent_var_senior_", patient, ".Rdata")
  load(deg.file)
  deg2<-deg.condition
  deg<-list(deg1, deg2)
  rm(deg.condition, deg1, deg2)
  # Group names
  groups<-c("Young", "Elderly", "MDS")
  res.file<-paste0("/home/mainciburu/scRNA/MDS_paper/GSEA/fgseares_", patient, ".txt")

  # Rearrange DEG tables
  CellTypes<-list()
  for(jj in 1:2){
    i<-sapply(deg[[jj]], is.data.frame)
    CellTypes[[jj]]<-deg[[jj]][i]
    names(CellTypes[[jj]])<-sapply(CellTypes[[jj]], function(x){unique(x$CellType)})
  }

  # GSEA
  fgseares<-data.frame()
  for(jj in 1:2){    # per comparison (Young vs MDS / Elderly vs MDS)
    fgseares.g<-data.frame()
    for(i in 1:length(CellTypes[[jj]])){   # per celltype
      print(i)
      celltype<-CellTypes[[jj]][[i]]
      rank<-celltype$avg_log2FC
      names(rank)<-rownames(celltype)
      fgseares.tmp<-fgsea::fgsea(pathways = m_list, stats = rank, nperm = 10000)
      fgseares.tmp$CellType<-names(CellTypes[[jj]])[i]
      fgseares.g<-rbind(fgseares.g, fgseares.tmp)
    }
    
    fgseares.g$GroupPair<-jj
    fgseares<-rbind(fgseares, fgseares.g)
  }

  fgseares<-fgseares[,c(1:3,5,9,10)]

  # Select significant results
  fgseares<-fgseares[fgseares$padj<=0.05,]

  # shorten pathway names
  rn<-fgseares$pathway
  rn<-sapply(rn, strsplit, split = "[_]")
  # x[-1] remove geneset source (hallmark, go...)
  rn<-sapply(rn, function(x){paste(x[-1], collapse = " ")})
  names(rn)<-NULL
  fgseares$pathway<-rn

  colnames(fgseares)<-c("GeneSet", "pval", "padj", "NES", "Contrast", "GroupPair")

  fgseares$Enriched.group[fgseares$GroupPair==1 & fgseares$NES>0]<-1
  fgseares$Enriched.group[fgseares$GroupPair==1 & fgseares$NES<0]<-3
  fgseares$Enriched.group[fgseares$GroupPair==2 & fgseares$NES>0]<-2
  fgseares$Enriched.group[fgseares$GroupPair==2 & fgseares$NES<0]<-3

  fgseares$Enriched.group<-factor(fgseares$Enriched.group)
  #fgseares$GeneSet<-factor(fgseares$GeneSet, levels = rev(fgseares$GeneSet[pathway_order]))
  #fgseares$Contrast<-factor(fgseares$Contrast, levels = rev(unique(fgseares$Contrast)))
  fgseares$Contrast<-factor(fgseares$Contrast, levels = unique(fgseares$Contrast))
  fgseares$GroupPair<-factor(fgseares$GroupPair, labels = c(1, 2), levels = c(1, 2))

  fgseares$Patient<-patient

  write.table(fgseares, file = res.file, quote = F, sep = "\t", row.names = F)
}

########################### Plot ###########################
res.path<-"/home/mainciburu/scRNA/MDS_paper/GSEA/"
mds1<-read.table(paste0(res.path, "fgseares_MDS1.txt"), header = T, sep = "\t")
mds2<-read.table(paste0(res.path, "fgseares_MDS2.txt"), header = T, sep = "\t")
mds3<-read.table(paste0(res.path, "fgseares_MDS3.txt"), header = T, sep = "\t")
mds4<-read.table(paste0(res.path, "fgseares_MDS4.txt"), header = T, sep = "\t")
groups<-c("Young", "Elderly", "MDS")

cols<-col.condition[names(col.condition)%in%groups]
names(cols)<-c("1", "2", "3")

fgseares<-rbind(mds1, mds2, mds3, mds4)
fgseares$Enriched.group<-factor(fgseares$Enriched.group)
fgseares$Contrast<-factor(fgseares$Contrast, levels = c("HSC", "LMPP", "GMP", "GMP_Granulocytes", 
                                                        "Monocytes", "pDC", "CLP", "MEP", "Megakaryocytes", 
                                                        "Erythroid_early", "Erythroid_late", "Basophils"))
fgseares$Patient<-factor(fgseares$Patient, levels = c("MDS1", "MDS2", "MDS3", "MDS4"))
fgseares$Group<-paste0(fgseares$Contrast, "_", fgseares$Patient)
lev<-paste0(rep(levels(fgseares$Contrast), each = 4), "_", 
            rep(levels(fgseares$Patient), length(levels(fgseares$Contrast))))
fgseares$Group<-factor(fgseares$Group, levels = lev)


### Young vs MDS
fgseares.y<-fgseares[fgseares$GroupPair==1,]

# Remove genesets with just one significant result
i<-table(fgseares.y$GeneSet)
i<-names(i)[i==1]
i<-!fgseares.y$GeneSet%in%i
fgseares.y<-fgseares.y[i,]

# Order genesets
enrichment_diff<-sapply(X = unique(fgseares.y$GeneSet), 
           function(x){
             a<-sum(fgseares.y[fgseares.y$GeneSet==x,"Enriched.group"]==1) 
             b<-sum(fgseares.y[fgseares.y$GeneSet==x,"Enriched.group"]==3) 
             return(c(a, b))
          }
)
colnames(enrichment_diff)<-unique(fgseares.y$GeneSet)
rownames(enrichment_diff)<-groups[c(1,3)]
i<-enrichment_diff[1,]>enrichment_diff[2,]
i<-which(i)
ix1<-order(enrichment_diff[1,], decreasing = T)
ix1<-ix1[ix1%in%i]

i<-enrichment_diff[1,]<=enrichment_diff[2,]
i<-which(i)
ix2<-order(enrichment_diff[2,], decreasing = F)
ix2<-ix2[ix2%in%i]
ix<-c(ix1, ix2)
if(length(ix)==length(unique(fgseares.y$GeneSet))){
    fgseares.y$GeneSet<-factor(fgseares.y$GeneSet, levels = unique(fgseares.y$GeneSet)[ix])
}

p<-ggplot(fgseares.y, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p<-p + theme_bw() + scale_color_manual(values = cols[c(1,3)], labels = groups[c(1,3)]) +
  scale_fill_manual(values = cols) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,6)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 16)) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched group") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
p<-p + theme(text = element_text(face = "bold"))

pdf("/home/mainciburu/scRNA/figures/supp_figure6/GSEA_hallmarks_latent_vars_young_MDS.pdf", width = 15, height = 14, useDingbats = F)
p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical") + ggtitle("Hallmark")
dev.off()

### Elderly vs MDS
fgseares.s<-fgseares[fgseares$GroupPair==2,]

# Remove genesets with just one significant result
i<-table(fgseares.s$GeneSet)
i<-names(i)[i==1]
i<-!fgseares.s$GeneSet%in%i
fgseares.s<-fgseares.s[i,]

# Order genesets
enrichment_diff<-sapply(X = unique(fgseares.s$GeneSet), 
           function(x){
             a<-sum(fgseares.s[fgseares.s$GeneSet==x,"Enriched.group"]==2) 
             b<-sum(fgseares.s[fgseares.s$GeneSet==x,"Enriched.group"]==3) 
             return(c(a, b))
          }
)
colnames(enrichment_diff)<-unique(fgseares.s$GeneSet)
rownames(enrichment_diff)<-groups[c(2,3)]
i<-enrichment_diff[1,]>enrichment_diff[2,]
i<-which(i)
ix1<-order(enrichment_diff[1,], decreasing = T)
ix1<-ix1[ix1%in%i]

i<-enrichment_diff[1,]<=enrichment_diff[2,]
i<-which(i)
ix2<-order(enrichment_diff[2,], decreasing = F)
ix2<-ix2[ix2%in%i]
ix<-c(ix1, ix2)
if(length(ix)==length(unique(fgseares.s$GeneSet))){
    fgseares.s$GeneSet<-factor(fgseares.s$GeneSet, levels = unique(fgseares.s$GeneSet)[ix])
}

p<-ggplot(fgseares.s, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p<-p + theme_bw() + scale_color_manual(values = cols[c(2,3)], labels = groups[c(2,3)]) +
  scale_fill_manual(values = cols) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,6)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 16)) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched group") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
p<-p + theme(text = element_text(face = "bold"))

pdf("/home/mainciburu/scRNA/figures/figure4/GSEA_hallmarks_latent_vars_elderly_MDS.pdf", width = 15, height = 14, useDingbats = F)
p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical") + ggtitle("Hallmark")
dev.off()

###################################
