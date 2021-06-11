
####################################################################
# fit gam
# imput
  # x - pseudotime 
  # y - expression of 1 gene
  # weights - cell probabilities for 1 branch
  # pred_x - pseudotime in 500 bins
# output
  # gene expression prediction in 500 bins
  # standar deviation
####################################################################

gam_fit_predict<-function(x, y, weights=NULL, pred_x=NULL){
  if(is.null(weights)){
    weights<-rep(1, length(x))
  }
  # select weights > 0
  use_inds = which(weights > 0)
  # fit model
  df<-data.frame(x = x[use_inds], y = y[use_inds])
  model = gam(formula = y~s(x), data=df,
                   weights=weights[use_inds])
  # Predict expression for 500 bins
  if(is.null(pred_x)){
    pred_x = x
  }
  y_pred = predict(object=model, newdata=data.frame(x=pred_x))
  # Standard deviations
  p = predict(model, newdata=data.frame(x=x[use_inds]))
  n<-length(use_inds)
  sigma = sqrt(sum((y[use_inds] - p) ** 2) / (n - 2))
  stds = sqrt(1 + 1 / n + (pred_x - mean(x)) ** 2 /
                   sum((x - mean(x)) ** 2)) * sigma / 2
  return(cbind(y_pred, stds))
}

##################################################
# Compute gene trends
# Input
  # branch_prob - cells x branches
  # Pseudotime - cells x pst value
  # gene_exprs - cells x selected genes
  # lineages - selected branches
  # ncores 
  # res_path - path to write gene trends files
# Output
  # results - list of lists (branches) with trends and std deviation
  # csv - one per branch and trends/std with selected genes x 500 bins
##################################################

compute_gene_trends<-function(branch_prob, pseudotime, 
                              gene_exprs, lineages, ncores=2, res_path)
{
  library(foreach)
  library(gam)
  library(doParallel)

  if(is.null(lineages)){   # Select every lineage
    lineages = colnames(branch_prob)
  }
  genes<-colnames(gene_exprs)
  results<-vector("list", length = length(lineages))
  names(results)<-lineages
  for(x in 1:length(results)){
    results[[x]]<-list()
    results[[x]][[1]]<-matrix(0, nrow = length(genes), ncol = 500)
    results[[x]][[2]]<-matrix(0, nrow = length(genes), ncol = 500)
    names(results[[x]])<-c("trends", "std")
  }
  
  for(branch in lineages){
    br_cells = rownames(branch_prob)[branch_prob[,branch]>0.7]
    bins = seq(from = 0, to = max(pseudotime[br_cells,1]), length.out = 500)
    colnames(results[[branch]][["trends"]])<-bins
    rownames(results[[branch]][["trends"]])<-genes
    colnames(results[[branch]][["std"]])<-bins
    rownames(results[[branch]][["std"]])<-genes
  }
  
  for(branch in lineages){
    weights<-branch_prob[,branch]
    bins<-colnames(results[[branch]][["trends"]])
    bins<-as.numeric(bins)
    registerDoParallel(cores=ncores)
    res<-foreach::foreach(i = 1:length(genes), .combine = "cbind")%dopar%{
      #print(i)
      rx<-gam_fit_predict(x = pseudotime[,1],
                      y = gene_exprs[,i],
                      weights = weights,
                      pred_x = bins)
      return(rx)
    }
    for(i in 1:length(genes)){
      r<-i*2
      results[[branch]][["trends"]][i,]<-res[,(r-1)]
      results[[branch]][["std"]][i,]<-res[,r]
    }
  }
  # write to feather
  for(branch in lineages){
    write.csv(x = results[[branch]][["trends"]], file = paste0(res_path, branch, "_trends.csv"))
    write.csv(x = results[[branch]][["std"]], file = paste0(res_path, branch, "_std.csv"))
  }
  return(results)
}


##################################################
# Read gene trend files
# Input
  # lineages - selected branches
  # res_path - path to write gene trends files
# Output
  # results - list of lists (branches) with trends and std deviation
##################################################

read_trends<-function(lineages, res.path){
  results<-vector(mode = "list", length = length(lineages))
  names(results)<-lineages
  for(branch in lineages){
    trend_file<-read.csv(paste0(res.path, branch, "_trends.csv"), row.names = 1)
    colnames(trend_file)<-gsub(pattern = "X", replacement = "", x = colnames(trend_file))
    std_file<-read.csv(paste0(res.path, branch, "_std.csv"), row.names = 1)
    colnames(std_file)<-colnames(trend_file)
    results[[branch]]<-vector(mode = "list", length = 2)
    names(results[[branch]])<-c("trends", "std")
    results[[branch]][["trends"]]<-trend_file
    results[[branch]][["std"]]<-std_file
  }
  return(results)
}
