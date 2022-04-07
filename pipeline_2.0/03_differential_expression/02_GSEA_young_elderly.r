######################
# GSEA per cluster between young and elderly
# MsigDb Hallmark pathways
# Plot results
##########################

library(fgsea)
library(ggplot2)
library(msigdbr)
library(foreach)
library(doParallel)
source("/home/mainciburu/scRNA/colors.r")


###### Input data ################
# differential expression results
deg.file<-"/home/mainciburu/scRNA/btwn_condition/deg_MAST_latent_var_young_senior.Rdata"
load(deg.file)

# Group names
groups<-c("Young", "Elderly")

# Hallmark genesets
m_df = msigdbr(species = "Homo sapiens")
m_df = m_df[m_df$gs_cat == "H",]
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
###############################

###### GSEA ###############

i<-sapply(deg.condition, is.data.frame)
CellTypes<-deg.condition[i]
names(CellTypes)<-sapply(CellTypes, function(x){unique(x$CellType)})

# Remove ProB, only present in elderly3
CellTypes<-CellTypes[-c(which(names(CellTypes)=="ProB"))]

fgseares<-data.frame()
for(i in 1:length(CellTypes)){
  print(i)
  celltype<-CellTypes[[i]]
  rank<-celltype$avg_log2FC
  names(rank)<-rownames(celltype)
  fgseares.tmp<-fgsea::fgsea(pathways = m_list, stats = rank, nperm = 10000)
  fgseares.tmp$CellType<-names(CellTypes)[i]
  fgseares<-rbind(fgseares, fgseares.tmp)
}

fgseares<-fgseares[,c(1:3,5,9)]

# Select significant results
fgseares<-fgseares[fgseares$padj<=0.05,]

# shorten pathway names
rn<-fgseares$pathway
rn<-sapply(rn, strsplit, split = "[_]")
rn<-sapply(rn, function(x){paste(x[-1], collapse = " ")})
names(rn)<-NULL
fgseares$pathway<-rn

colnames(fgseares)<-c("GeneSet", "pval", "padj", "NES", "Contrast")

fgseares$Enriched.group[fgseares$NES>0]<-1
fgseares$Enriched.group[fgseares$NES<0]<-2

write.table(fgseares, file = "/home/mainciburu/scRNA/btwn_condition/fgseares_MAST_latent_var_young_senior.txt", quote = F, sep = "\t", row.names = F)


###### Plot ########
fgseares<-read.table("/home/mainciburu/scRNA/btwn_condition/fgseares_MAST_latent_var_young_senior.txt", header = T, sep = "\t", stringsAsFactors=F)
fgseares$Enriched.group<-factor(fgseares$Enriched.group, levels = c(1, 2))
fgseares$Contrast<-factor(fgseares$Contrast, levels = names(CellTypes))

# Remove genesets with just one significant result
i<-table(fgseares$GeneSet)
i<-names(i)[i==1]
i<-!fgseares$GeneSet%in%i
fgseares<-fgseares[i,]

# Order genesets
enrichment_diff<-sapply(X = unique(fgseares$GeneSet), 
           function(x){
             a<-sum(fgseares[fgseares$GeneSet==x,"Enriched.group"]==1) 
             b<-sum(fgseares[fgseares$GeneSet==x,"Enriched.group"]==2) 
             return(c(a, b))
          }
)
colnames(enrichment_diff)<-unique(fgseares$GeneSet)
rownames(enrichment_diff)<-groups
i<-enrichment_diff[1,]>enrichment_diff[2,]
i<-which(i)
ix1<-order(enrichment_diff[1,], decreasing = T)
ix1<-ix1[ix1%in%i]

i<-enrichment_diff[1,]<enrichment_diff[2,]
i<-which(i)
ix2<-order(enrichment_diff[2,], decreasing = F)
ix2<-ix2[ix2%in%i]
ix<-c(ix1, ix2)
if(length(ix)==length(unique(fgseares$GeneSet))){
    fgseares$GeneSet<-factor(fgseares$GeneSet, levels = unique(fgseares$GeneSet)[ix])
}else{
    fgseares$GeneSet<-factor(fgseares$GeneSet, levels = c(unique(fgseares$GeneSet)[ix], unique(fgseares$GeneSet)[-ix]))
}

groups<-c("Young", "Elderly")
cols<-col.condition[names(col.condition)%in%groups]
names(cols)<-NULL

p<-ggplot(fgseares, aes(Contrast, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES))) +
  scale_shape_identity() 

p<-p + theme_bw() + scale_color_manual(values = cols, labels = groups) +
  scale_fill_manual(values = cols) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(2,8)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 16)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 16)) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched group") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
p<-p + theme(text = element_text(face = "bold"))

pdf("/home/mainciburu/scRNA/figures/figure1/gsea_hallmark_MAST_latent_var.pdf",
    width = 11.5, height = 10, useDingbats = F)
p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical") + ggtitle("Hallmark")
dev.off()
