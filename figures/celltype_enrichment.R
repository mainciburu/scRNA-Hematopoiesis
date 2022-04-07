library(enrichR)
library(ggplot2)
library(Seurat)
library(reshape2)
library(xlsx)
library(msigdbr)
library(tidyr)
library(dplyr)
library(plyr)
library(ggsignif)

########
# Tried with bioGPS and other .gmt from msigdb => no good results
########

######### custom genesets ##################
# Load cell type signatures
velten.sig<-read.xlsx("scRNA/btwn_condition/celltype_signatures_velten_2017.xlsx", 
                      sheetIndex = 1)
velten.sig$Gene<-as.character(velten.sig$Gene)
velten.sig$Gene<-substr(velten.sig$Gene, start = 1, stop = nchar(velten.sig$Gene) - 18)
velten.sig$ID<-1:nrow(velten.sig)
velten.sig<-velten.sig %>% spread(key = CellType, value = Gene) %>% select(-ID) %>% as.list()
velten.sig<-lapply(velten.sig, function(x) x[!is.na(x)])
names(velten.sig)<-paste0(names(velten.sig), "_velten")

kara.sig<-read.xlsx("scRNA/btwn_condition/lymphomyeloid_signatures_karamitros_2018.xlsx", 
                    sheetIndex = 1)
kara.sig$Gene<-as.character(kara.sig$Gene)
kara.sig$ID<-1:nrow(kara.sig)
kara.sig<-kara.sig %>% spread(key = CellType, value = Gene) %>% select(-ID) %>% as.list()
kara.sig<-lapply(kara.sig, function(x) x[!is.na(x)])
names(kara.sig)<-paste0(names(kara.sig), "_kara")

load("scRNA/markersBulk_progenitors.Rdata")
bulk.sig<-markersBulk_progenitors
bulk.sig<-bulk.sig[,c(1,4)]
bulk.sig<-bulk.sig[!is.na(bulk.sig$Gene),]
bulk.sig<-bulk.sig[bulk.sig$Gene!="",]
bulk.sig$ID<-1:nrow(bulk.sig)
bulk.sig<-bulk.sig %>% spread(key = CellType, value = Gene) %>% select(-ID) %>% as.list()
bulk.sig<-lapply(bulk.sig, function(x) x[!is.na(x)])
bulk.sig<-bulk.sig[1:5]
names(bulk.sig)<-paste0(names(bulk.sig), "_bulk")


list.sig<-c(velten.sig, kara.sig, bulk.sig)
term2gene<-data.frame()
for(i in 1:length(list.sig)){
  tmp<-cbind(names(list.sig)[i], list.sig[[i]])
  term2gene<-rbind(term2gene, tmp)
}
colnames(term2gene)<-c("TERM", "GENE")
#################################################

# Load cluster markers
deg.young<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_young.rds")
deg.senior<-readRDS("scRNA/btwn_celltype/deg_mast_celltype_senior.rds")
deg.young$Condition<-"Young"
deg.senior$Condition<-"Elderly"

# Load background
bg<-read.table("scRNA/common_genes.txt")
bg<-bg$V1

# Enrichments
res.young<-vector(mode = "list", length = length(levels(deg.young$cluster)))
res.senior<-vector(mode = "list", length = length(levels(deg.young$cluster)))

names(res.young)<-names(res.senior)<-levels(deg.young$cluster)
for(ii in 1:length(names(res.young))){
  cell<-names(res.young)[ii]
  print(paste0("-------", cell, "----------"))
  cl.markers.young<-deg.young[deg.young$cluster%in%cell & deg.young$p_val_adj<0.001,]
  cl.markers.senior<-deg.senior[deg.senior$cluster%in%cell & deg.senior$p_val_adj<0.001,]
  if(nrow(cl.markers.young)>0){
    res.tmp.young<-clusterProfiler::enricher(gene = cl.markers.young[,"gene"], pvalueCutoff = 0.1, 
                                             universe = bg, TERM2GENE = term2gene)
    res.young[[ii]]<-cbind(res.tmp.young@result, cell)
  }else{res.young[[ii]]<-NULL}
  if(nrow(cl.markers.senior)>0){
    res.tmp.senior<-clusterProfiler::enricher(gene = cl.markers.senior[,"gene"], pvalueCutoff = 0.1, 
                                             universe = bg, TERM2GENE = term2gene)
    res.senior[[ii]]<-cbind(res.tmp.senior@result, cell)
  }else{res.senior[[ii]]<-NULL}
}
res.young<-do.call("rbind", res.young)
res.senior<-do.call("rbind", res.senior)

colnames(res)[10]<-"CellType"

res.young<-res.young[res.young$p.adjust<0.1,]
res.senior<-res.senior[res.senior$p.adjust<0.1,]

res.young$Condition<-"Young"
res.senior$Condition<-"Senior"
res<-rbind(res.young, res.senior)
colnames(res)[10]<-"CellType"

# Plot
df<-data.frame(CellType=res$CellType, Condition=res$Condition, GeneSet=res$ID, 
               PValue=res$p.adjust, OddsRatio=res$GeneRatio)
or<-strsplit(as.character(df$OddsRatio), "/")
or<-lapply(or, as.numeric)
df$OddsRatio<-sapply(or, function(x){x[1]/x[2]})
df$CellType<-factor(df$CellType, levels = levels(res$CellType))
df$Condition<-factor(df$Condition, levels = c("Young", "Senior"))
df$GeneSet<-factor(df$GeneSet, levels = rev(c("HSC_bulk", "GMP_kara", "CMP_bulk", "CMP & GMP_bulk",
                                          "NeutP_velten", "EBM_velten", "MonoDC_velten",  
                                          "LMPP_kara", "CLP_bulk", "B_velten", "MEP_bulk", "MEP_velten")))
pdf("scRNA/figures_dec20/cluster_markers_enrichments.pdf", useDingbats = F, width = 10, height=7)
ggplot(df, aes(CellType, GeneSet, colour=-log10(PValue), size=OddsRatio)) + geom_point() + 
  facet_grid(~Condition) +
  theme_bw() + scale_color_gradient(low = "#FEE090", high = "#A50026") +
  theme(axis.text.y = element_text(hjust = 1, size = 16)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 16)) 
dev.off()
