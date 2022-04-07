#############################################
#########        Functions         ##########
#############################################


#######################################
#########   1) Integration   ##########
#######################################


# init_seurat - create seurat object and calculate % of mitochondrial genes
# input
	# name: sample name
	# path: path to the filtered 10X matrix
	# cols: string to add to the barcode names
	# condition: Young, Senior or MDS
# output: seurat object 

init_seurat<-function(name, path, cols, condition){
	init.data<-Read10X(data.dir = path)
	colnames(init.data)<-paste0(colnames(init.data), "_", cols)
	seurat.obj<-CreateSeuratObject(counts = init.data, min.cells = 3, min.features = 100, project = name)
	seurat.obj$Patient<-name
	seurat.obj$Condition<-condition

	mito.genes<-grep("MT-", rownames(seurat.obj), value = T)
	percent.mito<-Matrix::colSums(GetAssayData(seurat.obj, "counts")[mito.genes,])/Matrix::colSums(seurat.obj)
	seurat.obj$percent.mito<-percent.mito
	return(seurat.obj)
}

# preprocess_seurat - normalize, calculate cell cycle scores, scale, regress effects, find variable genes and PCA
# input
	# seurat.obj: object
	# npcs: number of principal component to calculate
	# vars: variables from metadata to regress
# output: processed seurat object

preprocess_seurat<-function(seurat.obj, npcs, vars){
	seurat.obj<-NormalizeData(seurat.obj)
	seurat.obj<-CellCycleScoring(object = seurat.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
	seurat.obj<-ScaleData(object = seurat.obj, vars.to.regress = vars)
	seurat.obj<-FindVariableFeatures(seurat.obj, selection.method = "vst",
                            nfeatures = 2000)
	seurat.obj<-RunPCA(seurat.obj, npcs = npcs)
	return(seurat.obj)
}

# cluster_seurat: cluster cells with specific resolution, plot and calculate average silhouette
# input
	# seurat.obj: object
	# npcs: number of significant components
	# resolution: clustering resolution
# output: seurat object with clusters, cluster umap plot and silhouette barplot

cluster_seurat<-function(seurat.obj, npcs, resolution){
	seurat.obj<-FindNeighbors(object = seurat.obj, reduction = "pca", dims = 1:npcs)
	seurat.obj<-FindClusters(seurat.obj, resolution = resolution)
	col<-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.obj$seurat_clusters)))
	print(DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", cols = col))
	# Cluster silhouette
	d<-dist(x = seurat.obj@reductions$pca@cell.embeddings[,1:npcs])
	s<-silhouette(x = as.numeric(seurat.obj$seurat_clusters), dist = d)
	summary(s)
	s.avg<-as.numeric(summary(s)$clus.avg.widths)
	c<-length(unique(seurat.obj$seurat_clusters)) - 1
	barplot(s.avg, horiz = T, names.arg = as.character(0:c), col = col)
	print("Silhouette per cluster")
	names(s.avg)<-as.character(0:c)
	print(s.avg)
	print(paste0("Average silhouette: ", mean(s.avg)))
	return(seurat.obj)
}

# vlnplot- create umap + vlnplot with lineage markers per cluster
# input
	# seurat.obj: object
	# title: plot title
	# clusters: meta.data column to make groups
# output: plot in pdf

vlnplot<-function(seurat.obj, cluster, title=NULL, col=NULL){
  InfoData<-data.frame(x=seurat.obj@reductions$umap@cell.embeddings[,1], 
                       y=seurat.obj@reductions$umap@cell.embeddings[,2],
                       Cluster=seurat.obj@meta.data[,cluster])
  # Markers
  hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
  lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
  cc<-c("CDC20", "TOP2A")
  gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
  granul<-c("ELANE", "AZU1", "CEBPA", "CEBPE", "CST7")
  mono<-c("LYZ", "CSTA")
  dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4")
  t<-c("JCHAIN", "IKZF1", "CYTH1")
	nk<-c("TSC22D1", "CXXC5", "HOXA9", "HOXA10")
  clp<-c("IL7R", "DNTT")
  prob<-c("VPREB1", "EBF1", "CD79A", "CD79B")
  mep<-c("NFE2", "HFS1", "TAL1")
  mk<-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
  ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
  baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")
  
  markers<-c(hsc, lmpp, cc, gmp, granul, mono, dc, clp, prob, t, mep, mk, ery, baso)
  markers<-markers[markers%in%rownames(seurat.obj)]
  
  InfoData<-cbind(InfoData,t(as.matrix(seurat.obj@assays$RNA@data[markers,])))
  InfoData2<-melt(data = InfoData, measure.vars = markers)
  if(is.null(col)){
    col<-colorRampPalette(brewer.pal(12, "Paired"))(nrow(unique(seurat.obj[[cluster]])))
  }
  Pvln<-ggplot(InfoData2, aes(x=Cluster, y=value, fill = Cluster))+facet_grid(variable~.)+geom_violin(scale = "width")+ theme(legend.position="none") + labs(x=NULL, y = "Counts")+scale_fill_manual(values=col) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(strip.text.y = element_text(angle = 0)) + theme(strip.background = element_blank()) + scale_y_continuous(limits = c(0,6))
  Pumap<-DimPlot(seurat.obj, reduction = "umap", group.by = cluster, pt.size = 0.5, cols = col) + labs(colour = "Cluster") + ggtitle(label = title)
  plot_grid(Pumap, Pvln, ncol = 1, rel_heights = c(1, 1.5))
}


#######################################
#########   4) Trajectories   #########
#######################################

plot_trend_2conditions <- function(data, gene){
  test_results <- wilcox.test(data$trend ~ data$condition)
  tv_value <- sum(abs(data$trend[data$condition=="young"] - data$trend[data$condition=="senior"]))
    p <- ggplot(data, aes(x = time, y = trend, group = condition, color=condition)) + 
       geom_ribbon(aes(ymin = ymin , ymax = ymax), linetype = 2, alpha=0.1)+
       geom_line(size = 1.5) + 
       scale_color_manual(values = c('young' = '#C9453F', 'senior'='#557CAE')) +
       labs(title = gene,
            subtitle = paste0('Wilcox p-value:', signif(test_results$p.val, 3),
                                '; Total Variation:', signif(tv_value,3)),
            x = "Pseudotime", y = "Expression")
    p<-p + theme(text = element_text(face = "bold", color="black"), 
                           plot.title = element_text(size = 28), 
                           plot.subtitle = element_text(size = 16), 
                           legend.position = "none",
                           axis.text = element_text(size = 30, color="black"),
                           axis.title = element_text(size = 20),
                           axis.line = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.35, "cm")) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())
    return(p)
}


plot_trend_multiple_conditions <- function(data, gene){
  #test_results <- wilcox.test(data$trend ~ data$condition)
  #tv_value <- sum(abs(data$trend[data$condition=="young"] - data$trend[data$condition=="senior"]))
    p <- ggplot(data, aes(x = time, y = trend, group = condition, color=condition)) + 
       geom_ribbon(aes(ymin = ymin , ymax = ymax), linetype = 2, alpha=0.1)+
       geom_line(size = 1.5) + 
       scale_color_manual(values = col.condition) +
       labs(title = gene,
            #subtitle = paste0('Wilcox p-value:', signif(test_results$p.val, 3),
            #                    '; Total Variation:', signif(tv_value,3)),
            x = "Pseudotime", y = "Expression")
    p<-p + theme(text = element_text(face = "bold", color="black"), 
                           plot.title = element_text(size = 28), 
                           plot.subtitle = element_text(size = 16), 
                           #legend.position = "none",
                           axis.text = element_text(size = 30, color="black"),
                           axis.title = element_text(size = 20),
                           axis.line = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.35, "cm")) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())
    return(p)
}


add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                   grobs = new.flag,
                                   t = 4, 
                                   l = 4
  )

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid.newpage()
  grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}



#######################################
###########     5) GRN     ############
#######################################

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
