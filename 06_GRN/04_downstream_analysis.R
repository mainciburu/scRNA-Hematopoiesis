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


BinarizeAUC<-function(aucMatrix,thrP = 0.01,smallestPopPercent = 0.25,plotHist = TRUE,densAdjust = 2,nBreaks = 100)
{
  if (any(rowSums(aucMatrix) == 0)) 
  {
    warning("Skipping genesets with all AUC 0: ", paste(rownames(aucMatrix)[which(rowSums(aucMatrix) == 0)], collapse = ", "), immediate. = TRUE)
    aucMatrix <- aucMatrix[rowSums(aucMatrix) > 0, , drop = FALSE]
  }
  
  
  gSetName <- NULL
  
  assigment <- lapply(rownames(aucMatrix), function(gSetName) {
    print(gSetName)
    aucThr <- AUCell:::.auc_assignmnetThreshold_v6(aucRow = aucMatrix[gSetName,, drop = FALSE], 
                                                   thrP = thrP, 
                                                   smallestPopPercent = smallestPopPercent,
                                                   plotHist = plotHist, 
                                                   densAdjust = densAdjust,
                                                   nBreaks = nBreaks)
    
    assignedCells <- NULL
    BinVec<-rep(0,ncol(aucMatrix))
    
    if (!is.null(aucThr))
    {
      auc <- aucMatrix[gSetName, ]
      assignedCells <- names(auc)[which(auc > aucThr$selected)]
      assignedCells<-match(assignedCells,colnames(aucMatrix))
      BinVec[assignedCells]<-1
      Result<-list(aucThr = aucThr$selected, assignment = BinVec)
    }else{
      Result<-list(aucThr = NULL, assignment = BinVec)
    }
    
  })
  
  names(assigment) <- rownames(aucMatrix)
  
  Thresholds <- sapply(assigment, function(x) x[[1]])
  names(Thresholds)<-rownames(aucMatrix)
  Cells<- lapply(assigment, function(x) x[[2]])
  BinMat<-do.call(rbind,Cells)
  rownames(BinMat)<-rownames(aucMatrix)
  colnames(BinMat)<-colnames(aucMatrix)
  
  Binarized<-list(Thresholds=Thresholds,BinaryMatrix=BinMat)
  return(Binarized)
  
}


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

# Read regulons identified with PyScenic
Regulons<-ReadPyScenicGMT(file = "/home/mainciburu/scRNA/scenic/regulons/AML_regulons.gmt")

# Read RSS Scores
RSS<-read.delim(file="/home/mainciburu/scRNA/scenic/regulons/AML_rssCelltype.csv",sep=",")

# Read MetaData 
MetaData<-read.delim(file="/home/mainciburu/scRNA/scenic/data/aml_CellType.txt",sep=" ",header = T)
MetaData$CellType<-factor(MetaData$CellType,levels = c("HSC", "LMPP", "GMP", "GMP_Granulocytes", 
                                                 "Monocytes", "pDC", "CLP", "T_NK", "ProB", 
                                                 "MEP", "Megakaryocytes", "Erythroid_early", 
                                                 "Erythroid_late", "Basophils", "not assigned"))

# Read Regulon AUC matrix
AUC<-read.csv(file="/home/mainciburu/scRNA/scenic/regulons/AML_regulons_AUCMat.csv",stringsAsFactors = F)
load("/home/mainciburu/scRNA/scenic/results/binAUC_AML.Rdata")

# We will only work with Activation Regulons
iix<-which(Regulons$Mode=="Activation")
Regulons<-Regulons[iix,]
Regulons[,1]<-gsub("[(+)]","",Regulons[,1])

# Only Keep Activation and Cells Included in the Seurat Object
Regs<-AUC[,1]
Regs<-gsub("[(+)]","",Regs)
rownames(AUC)<-Regs
AUC<-AUC[,2:ncol(AUC)]
AUC<-as.matrix(AUC)
AUC<-AUC[,rownames(MetaData)]

## Transform RSS Matrix
rownames(RSS)<-gsub("[(+)]","",RSS[,1])
RSS<-RSS[,2:ncol(RSS)]

#### Only Correct If Required
AUC<-removeBatchEffect(AUC,batch = MetaData$Donor)

## Binarize AUC
BinAUC<-BinarizeAUC(AUC,plotHist = FALSE)
BinaryAUC<-BinAUC$BinaryMatrix
save(BinaryAUC, file = "/home/mainciburu/scRNA/scenic/results/binAUC_AML.Rdata")
load("/home/mainciburu/scRNA/scenic/results/binAUC_AML.Rdata")

## Create Regulon Activity Per Cell Type (average per cluster)
Dx<-split(rownames(MetaData),MetaData$CellType)
regulonActivity_byCellType<-sapply(Dx[as.numeric(summary(Dx)[,1])>15],
                                   function(cells){rowSums(BinaryAUC[,cells])/length(cells)})
rownames(regulonActivity_byCellType)<-rownames(BinaryAUC)


## Select Top 5 Regulons per CellType
TFs<-lapply(RSS,function(X){return(rownames(RSS)[order(X,decreasing = T)][1:5])})
TFs<-TFs[levels(MetaData$CellType)]

for(jj in colnames(regulonActivity_byCellType))
{
  
  TFs[[jj]]<-names(sort(regulonActivity_byCellType[TFs[[jj]],jj],decreasing = T))
  
}

TFs<-TFs[levels(MetaData$CellType)]

# save TF list
x<-unique(do.call("c", TFs))
write.table(x, file = "/home/mainciburu/scRNA/scenic/results/AML_TFs.txt", sep = "\t", col.names = F)

regulonActivity_byCellType<-regulonActivity_byCellType[,levels(MetaData$CellType)[levels(MetaData$CellType)%in%colnames(regulonActivity_byCellType)]]
regulonActivity_byCellType<-regulonActivity_byCellType[,colnames(regulonActivity_byCellType)!="not assigned"]

Info<-data.frame(CellType=colnames(regulonActivity_byCellType))
rownames(Info)<-Info[,1]
Info$CellType<-factor(Info$CellType,levels = Info$CellType)

MyCols<-c(col.young.v3[names(col.young.v3)%in%Info$CellType], "not assigned" = "grey")

annotation_colors = list(CellType=MyCols)

pheatmap(regulonActivity_byCellType[unique(unlist(TFs[names(TFs)%in%colnames(regulonActivity_byCellType)])),],
         cluster_rows = F,cluster_cols = F, annotation_col=Info ,annotation_colors = annotation_colors,
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(75),border_color = NA,main = "AML",fontsize_row = 8,
         filename="/home/mainciburu/scRNA/scenic/plots/AML_heatmap.pdf")


############# Enrichments #########################
Rx<-as.data.frame(unlist(TFs[names(TFs)%in%colnames(regulonActivity_byCellType)]))
Rx<-data.frame(CellType=rownames(Rx),Regulon=Rx[,1])
Rx$CellType<-gsub("[0-9]","",Rx$CellType)
Rx$Target<-as.character(Regulons[match(Rx$Regulon,Regulons$Regulon),"Targets"])

DD<-strsplit(Rx$Target,",")
DataEnrichment<-data.frame(gene=unlist(DD),cluster=rep(Rx$CellType,sapply(DD, length)))

## Remove CLP, ProB, T_NK and Monocytes (<15 cells)
DataEnrichment<-DataEnrichment[!DataEnrichment$cluster%in%c("CLP", "ProB", "T_NK", "Monocytes"),]

Res<-bitr(DataEnrichment[,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DataEnrichment$Entrez<-Res[match(DataEnrichment$gene,Res[,1]),2]
DataEnrichment$Group<-"AML"

GO <- compareCluster(split(DataEnrichment$Entrez,DataEnrichment$cluster) ,fun="enrichGO",ont="BP", OrgDb='org.Hs.eg.db',pvalueCutoff=0.05)
GO@compareClusterResult$Cluster<-factor(GO@compareClusterResult$Cluster,levels=levels(MetaData$label))
GO@compareClusterResult<-GO@compareClusterResult[order(GO@compareClusterResult$Cluster,GO@compareClusterResult$p.adjust,decreasing = F),]
GO <- setReadable(GO, 'org.Hs.eg.db', 'ENTREZID')

DataPlot<-GO@compareClusterResult

DataPlot$GenRat<-sapply(strsplit(DataPlot$GeneRatio,"/"),function(X){as.numeric(X[1])/as.numeric(X[2])})

iix<-c()
for(jj in unique(DataPlot$Cluster))
{
  
  AA<-DataPlot[which(DataPlot$Cluster==jj),]
  Trms<-AA[1:5,"Description"]
  
  for(ii in 1:5)
  {
    iix<-c(iix,which(DataPlot$Description==Trms[ii]))
  }
}

GGData<-DataPlot[iix,]
GGData<-GGData[order(GGData$Cluster,GGData$p.adjust),]
GGData$Description<-factor(GGData$Description,levels=unique(GGData$Description))
index<-which(GGData$Cluster=="HSC"|GGData$Cluster=="LMPP"|GGData$Cluster=="GMP"|GGData$Cluster=="GMP_Granulocytes"|GGData$Cluster=="MEP")
pp<-ggplot(GGData[index,],aes(x=Description,y=-log10(pvalue),fill=Cluster))+
  geom_bar(stat="identity",position = "dodge")+coord_flip()+facet_grid(~Cluster)+xlab("-")+ylab("-log10 pvalue")+
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#006D2C","#B2DF8A","#FB9A99"))
pdf("scRNA/scenic/plots/AML_enrichment_barplot.pdf", width = 15, height = 6)
print(pp)
dev.off()


#### RSS For Supplamentary
CellOrder<-levels(MetaData$label) 
CellOrder<-CellOrder[!CellOrder%in%c("CLP", "ProB", "T_NK", "Monocytes", "not assigned")]

RSS<-read.delim(file="/home/mainciburu/scRNA/scenic/regulons/AML2_rssCelltype.csv",sep=",")
RSS[,1]<-gsub("[(+)]","",RSS[,1])
colnames(RSS)[1]<-"Regulon"
RSS<-melt(RSS)
RSS<-RSS[RSS$variable%in%CellOrder,]
RSS<-RSS %>% group_by(variable) %>% arrange(desc(value),.by_group = T)
colnames(RSS)<-c("Regulon","CellType","RSS")
RSS$CellType<-factor(RSS$CellType,levels = CellOrder)
NRegs<-length(unique(RSS$Regulon))
RSS$Position<-rep(1:NRegs,length(CellOrder))

Plot1<-ggplot(RSS[RSS$CellType%in%CellOrder[1:5],],aes(x=Position,y=RSS))+geom_point()+facet_grid(~CellType)+
  geom_label_repel(aes(label=ifelse(Position<=5,as.character(Regulon),'')),fill="lightgray",box.padding   = 0.35, 
                   point.padding = 0.5,max.overlaps = Inf,na.rm = T,size=2.5,fontface="bold")+theme_classic()+
  theme(strip.text.x = element_text(face="bold"))+ylim(c(0,0.7))

Plot2<-ggplot(RSS[RSS$CellType%in%CellOrder[6:10],],aes(x=Position,y=RSS))+geom_point()+facet_grid(~CellType)+
  geom_label_repel(aes(label=ifelse(Position<=5,as.character(Regulon),'')),fill="lightgray",box.padding   = 0.35, 
                   point.padding = 0.5,max.overlaps = Inf,na.rm = T,size=2.5,fontface="bold")+theme_classic()+
  theme(strip.text.x = element_text(face="bold"))+ylim(c(0,0.7))

pdf(file="scRNA/scenic/plots/AML_RSS.pdf",width=8.6,height=6, useDingbats = F)
Plot1+Plot2+plot_layout(ncol = 1)
dev.off()

