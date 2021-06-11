#### Files for GRN visualization in Cytoscape
library(Seurat)
library(igraph)
library(AUCell)
library(limma)
elibrary(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(ggrepel)
library(patchwork)

source("/home/mainciburu/scRNA/colors.r")

ReadPyScenicGMT<-function(file)
{
  Info<-readLines(file)
  
  MyResult<-vector("list",length=length(Info))
  for(jj in 1:length(Info))
  {
    A<-strsplit(Info[jj],",")
    Reg<-A[[1]][1]
    Mode<-ifelse(strsplit(Reg,"[()]")[[1]][2] == "+","Activation","Repression")
    Score<-as.numeric(gsub("score=","",A[[1]][3]))
    Targets<-paste(A[[1]][4:length(A[[1]])],collapse = ",")
    Res<-data.frame(Regulon=Reg,Mode=Mode,Score=Score,Targets=Targets)
    MyResult[[jj]]<-Res
    
  }
  
  Result<-do.call(rbind,MyResult)
  return(Result)
  
}
#################################################

# Read regulons identified with PyScenic
Regulons<-ReadPyScenicGMT(file = "/home/mainciburu/scRNA/scenic/regulons/AML_regulons.gmt")

# We will only work with Activation Regulons
iix<-which(Regulons$Mode=="Activation")
Regulons<-Regulons[iix,]
Regulons[,1]<-gsub("[(+)]","",Regulons[,1])

# Read RSS Scores
RSS<-read.delim(file="/home/mainciburu/scRNA/scenic/regulons/AML_rssCelltype.csv",sep=",")

## Transform RSS Matrix
rownames(RSS)<-gsub("[(+)]","",RSS[,1])
RSS<-RSS[,2:ncol(RSS)]

# Read Regulon AUC matrix
load("/home/mainciburu/scRNA/scenic/results/binAUC_AML.Rdata")
BinaryAUC<-BinAUC$BinaryMatrix

AllTFs<-intersect(rownames(BinaryAUC),rownames(RSS))

## Define Same Set of Regulons
BinaryAUC<-BinaryAUC[AllTFs,]
Regulons<-Regulons[match(AllTFs,Regulons[,1]),]

## Select Top 5 Regulons per CellType
TFs<-lapply(RSS,function(X){return(rownames(RSS)[order(X,decreasing = T)][1:5])})
TFs<-TFs[CellOrder]


# Read Adjacency and TF list
# Adjacency file
Adj<-read.delim(file="/home/mainciburu/scRNA/scenic/adjacencies/AML_adjacency.tsv",sep="\t",header = T)

# Load TFs
AllTFs<-read.delim(file="/home/mainciburu/scRNA/scenic/results/AML_TFs.txt",sep="\t",header = F, stringsAsFactors=F)

TFs<-data.frame(CellType=rep(na.omit(names(TFs)),each=5),TF=unlist(TFs))

Regulons<-Regulons[match(TFs$TF,Regulons$Regulon),]
Regulons$CellType<-TFs[match(Regulons$Regulon,TFs$TF),1]

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
Nodes[Nodes$NAME%in%AllTFs$V1,2]<-"TF"

Dupls<-names(which(table(Network$TARGET)>1))
Nodes$CellType<-as.character(Nodes$CellType)
Nodes$CellType[Nodes$NAME%in%Dupls]<-"Multiple"


# Write files for cytoscape
write.table(Network,file="/home/mainciburu/scRNA/scenic/results/NetworkAML.txt",
            sep="\t",col.names = T,row.names = F,quote = F)
write.table(Nodes,file="/home/mainciburu/scRNA/scenic/results/NodesAML.txt",
            sep="\t",col.names = T,row.names = F,quote = F)

