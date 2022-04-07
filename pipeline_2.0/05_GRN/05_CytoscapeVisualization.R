#### Files for GRN visualization in Cytoscape
#library(Seurat)
library(igraph)
library(AUCell)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(ggrepel)
library(patchwork)

source("/home/mainciburu/scRNA/colors.r")
source("/home/mainciburu/scRNA/pipeline_2.0/functions.r")

MySample<-"MDS1"
# Read regulons identified with PyScenic
Regulons<-ReadPyScenicGMT(file = paste0("/home/mainciburu/scRNA/scenic/regulons/", MySample, "_regulons.gmt"))

# We will only work with Activation Regulons
iix<-which(Regulons$Mode=="Activation")
Regulons<-Regulons[iix,]
Regulons[,1]<-gsub("[(+)]","",Regulons[,1])

# Read RSS Scores
RSS<-read.delim(file=paste0("/home/mainciburu/scRNA/scenic/regulons/", MySample, "_rssCelltype.csv"),sep=",")

## Transform RSS Matrix
rownames(RSS)<-gsub("[(+)]","",RSS[,1])
RSS<-RSS[,2:ncol(RSS)]

# Read Regulon AUC matrix
load(paste0("/home/mainciburu/scRNA/scenic/results/binAUC_", MySample, ".Rdata"))
#BinaryAUC<-BinAUC$BinaryMatrix

AllTFs<-intersect(rownames(BinaryAUC),rownames(RSS))

## Define Same Set of Regulons
BinaryAUC<-BinaryAUC[AllTFs,]
Regulons<-Regulons[match(AllTFs,Regulons[,1]),]

# Read Adjacency and TF list
# Adjacency file
Adj<-read.delim(file=paste0("/home/mainciburu/scRNA/scenic/adjacencies/", MySample, "_adjacency.tsv"),sep="\t",header = T)

# Load TFs
TFs<-read.delim(file=paste0("/home/mainciburu/scRNA/scenic/results/", MySample, "_TFs.txt"),sep="\t",header = F, stringsAsFactors=F)

colnames(TFs)<-c("TF", "CellType")

Regulons<-Regulons[match(TFs$TF,Regulons$Regulon),]
Regulons$CellType<-TFs[match(Regulons$Regulon,TFs$TF),2]

Targets<-strsplit(as.character(Regulons$Targets),",")
names(Targets)<-Regulons$Regulon

Network<-data.frame(SOURCE=rep(names(Targets),times=sapply(Targets,length)),
                    TARGET=unlist(Targets),
                    INTERACTION="Activation",
                    CellType=rep(Regulons$CellType,times=sapply(Targets,length)))


Network$IMPORTANCE<-Adj[match(paste0(Network$SOURCE,"-",Network$TARGET),paste0(Adj$TF,"-",Adj$target)),3]

Network<-Network[!is.na(Network$IMPORTANCE),]

Network<-unique(Network)

## Filter targets - remove those with importance < 3rd quartile
Idx<-c()
for(jj in unique(Network$SOURCE))
{
  
  AA<-Network[Network$SOURCE==jj,]
  Idx<-c(Idx,rownames(AA)[AA$IMPORTANCE>=quantile(AA$IMPORTANCE,0.75)])
}

Network<-Network[Idx,]

Nodes<-data.frame(NAME=unique(c(as.character(Network$SOURCE),as.character(Network$TARGET))),
				  TYPE="Target", stringsAsFactors=F)
Nodes[Nodes$NAME%in%Regulons$Regulon,2]<-"TF"


Nodes<-split(Nodes,Nodes$TYPE)
Nodes$TF$CellType<-Regulons[match(Nodes$TF$NAME,Regulons$Regulon),"CellType"]
Nodes$Target$CellType<-Network[match(Nodes$Target$NAME,Network$TARGET),"CellType"]
Nodes<-do.call(rbind,Nodes)
Nodes[Nodes$NAME%in%TFs$TF,2]<-"TF"

Dupls<-names(which(table(Network$TARGET)>1))
Nodes$CellType<-as.character(Nodes$CellType)
Nodes$CellType[Nodes$NAME%in%Dupls]<-"Multiple"


# Write files for cytoscape
write.table(Network,file=paste0("/home/mainciburu/scRNA/scenic/results/Network_", MySample, ".txt"),
            sep="\t",col.names = T,row.names = F,quote = F)
write.table(Nodes,file=paste0("/home/mainciburu/scRNA/scenic/results/Nodes_", MySample, ".txt"),
            sep="\t",col.names = T,row.names = F,quote = F)

