######################
# GSEA per cluster between conditions
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
deg.file<-"/home/mainciburu/scRNA/btwn_condition/deg_MAST_young_senior.Rdata"
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

fgseares<-data.frame()
for(i in 1:length(CellTypes)){
  print(i)
  celltype<-CellTypes[[i]]
  rank<-celltype$avg_logFC
  names(rank)<-rownames(celltype)
  fgseares.tmp<-fgsea::fgsea(pathways = m_list, stats = rank, nperm = 10000)
  fgseares.tmp$CellType<-names(CellTypes)[i]
  fgseares<-rbind(fgseares, fgseares.tmp)
}

# Adjust p.val to all comparisons
fgseares$padj2<-p.adjust(fgseares$pval, method = "BH")

fgseares<-fgseares[,c(1:3,5,9,10)]

# Select significant results
fgseares<-fgseares[fgseares$padj2<=0.05,]

# shorten pathway names
rn<-fgseares$pathway
rn<-sapply(rn, strsplit, split = "[_]")
rn<-sapply(rn, function(x){paste(x[-1], collapse = " ")})
names(rn)<-NULL
fgseares$pathway<-rn

colnames(fgseares)<-c("GeneSet", "pval", "padj", "NES", "Contrast", "padj2" )

fgseares$Enriched.group[fgseares$NES>0]<-1
fgseares$Enriched.group[fgseares$NES<0]<-2

fgseares$Enriched.group<-factor(fgseares$Enriched.group)
fgseares$GeneSet<-factor(fgseares$GeneSet)
#fgseares$Contrast<-factor(fgseares$Contrast, levels = rev(unique(fgseares$Contrast)))
fgseares$Contrast<-factor(fgseares$Contrast, levels = unique(fgseares$Contrast))

write.table(fgseares, file = "/home/mainciburu/scRNA/btwn_condition/fgseares_MAST.txt", quote = F, sep = "\t", row.names = F)
