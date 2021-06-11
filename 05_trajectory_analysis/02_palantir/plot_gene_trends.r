#########################################################

# r1 -> results set
# r2=NULL -> results 2 set
# gg -> genes to plot
# bb -> branches to plot 
# Modes (condition, gene, trajectory)
  # one condition, one genes, multiple trajectory -> 112
  # one condition, multiple genes, one trajectory -> 121
  # two conditions, one genes, one trajectory -> 211
  # three conditions, one genes, one trajectory -> 311
# col1 -> color set
# col2=NULL
# c1 -> condition 1 
# c2 -> condition 2

plot_trends<-function(r1, r2=NULL, r3=NULL, gg, bb=NULL, mode, 
                      col1=NULL, col2=NULL, col3=NULL, c1=NULL, c2=NULL, c3=NULL){

##### One condition, one gene, multiple trajectories --------------------
if(mode==112){
  
  p<-list()
  r<-1
  results<-r1
  col<-col1
  for(g in gg){
    df<-data.frame()
    i<-1
    for(branch in results){
      dfb<-data.frame(x=as.numeric(colnames(branch$trends)),
                      exprs=as.numeric(branch$trends[g,]),
                      branch=rep(names(results)[i], ncol(branch$trends)),
                      std=as.numeric(branch$std[g,]))
      df<-rbind(df, dfb)
      i<-i+1
    }
    #df$exprs[df$exprs<0]<-0
    df$ymin<-df$exprs - df$std
    df$ymax<-df$exprs + df$std
    p[[r]]<-ggplot(df, aes(x, exprs, group = branch, colour = branch)) + 
      geom_line() + theme(axis.text.x = element_blank()) +
      geom_ribbon(aes(x = x, ymin=ymin, ymax=ymax, fill = branch), linetype=0, alpha=0.1) +
      scale_color_manual(values = col) +     scale_fill_manual(values = col) + 
      labs(x = "Pseudotime", y = "Expression") + ggtitle(g) + theme(legend.position = "none")
    p[[r]]<-p[[r]] + theme(text = element_text(face = "bold"), title = element_text(size = 20), axis.text = element_text(size = 16)) +
                     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"))
    r<-r+1
  }
  return(plot_grid(plotlist = p, nrow = 1))
}
###### One condition, multiple genes, one trajectory --------------------

if(mode==121){
  results<-r1
  p<-list()
  r<-1
  for(b in bb){
    t<-results[[b]]$trends
    t.std<-results[[b]]$std
    
    df<-data.frame(x = as.numeric(colnames(t)))
    df<-cbind(df, t(t[gg,]))
    df<-melt(df, measure.vars = gg)
    df$std<-0
    std<-c()
    for(g in gg){std<-c(std, as.numeric(t.std[g,]))}
    df$std<-std
    df$ymin<-df$value - df$std
    df$ymax<-df$value + df$std
    df$ymin[df$ymin<0]<-0
    p[[r]]<-ggplot(df, aes(x = x, y = value, group = variable, colour = variable)) + 
      geom_line() + theme(axis.text.x = element_blank()) + 
      geom_ribbon(aes(x = x, ymin=ymin, ymax=ymax, fill = variable), linetype=0, alpha=0.1) + 
      labs(x = "Pseudotime", y = "Expression") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    r<-r+1
  }
  return(plot_grid(plotlist = p))
}

##### Two conditions, one gene, one trajectory --------------------------

if(mode==211){
  t1<-r1[[bb]]$trends
  t2<-r2[[bb]]$trends
  t1.std<-r1[[bb]]$std
  t2.std<-r2[[bb]]$std
  col<-c(col1, col2)
  names(col)<-c(c1, c2)
  p<-list()
  r<-1
  for(g in gg){
    df<-data.frame(x = as.numeric(c(colnames(t1), colnames(t2))),
                   Condition = c(rep(c1, ncol(t1)), rep(c2, ncol(t2))))
    df$Condition<-factor(df$Condition, levels=c(c1, c2))
    df$exprs<-c(as.numeric(t1[g,]), as.numeric(t2[g,]))
    df$std<-c(as.numeric(t1.std[g,]), as.numeric(t2.std[g,]))
    df$ymin<-df$exprs - df$std
    df$ymax<-df$exprs + df$std
    
    p[[r]]<-ggplot(df, aes(x, exprs, group = Condition, colour = Condition)) +
      geom_line(size = 1.5) + theme(axis.text.x = element_blank()) + 
      scale_color_manual(values = col) +
      geom_ribbon(aes(x = x, ymin=ymin, ymax=ymax), linetype="dotted", alpha=0.1) + 
      labs(x = "Pseudotime", y = "Expression") + 
      ggtitle(g)
    p[[r]]<-p[[r]] + theme(text = element_text(face = "plain", color="black"), 
                           title = element_text(size = 36), 
                           legend.text = element_text(size = 32),
                           axis.text.y = element_text(size = 30, color="black"),
                           axis.text.x = element_text(size=30, angle = 30, hjust = 1, color="black"),
                           axis.line = element_line(colour = "black", size = 1.5),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.35, "cm")) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())
    r<-r+1
  }
  if(r==2){return(p[[1]])}
  if(r>2){return(plot_grid(plotlist = p))}
}


##### THREE conditions, one gene, one trajectory --------------------------

if(mode==311){
  t1<-r1[[bb]]$trends
  t2<-r2[[bb]]$trends
  t3<-r3[[bb]]$trends
  t1.std<-r1[[bb]]$std
  t2.std<-r2[[bb]]$std
  t3.std<-r3[[bb]]$std
  col<-c(col1, col2, col3)
  names(col)<-c(c1, c2, c3)
  p<-list()
  r<-1
  for(g in gg){
    df<-data.frame(x = as.numeric(c(colnames(t1), colnames(t2), colnames(t3))),
                   Condition = c(rep(c1, ncol(t1)), rep(c2, ncol(t2)),
                                 rep(c3, ncol(t3))))
    df$Condition<-factor(df$Condition, levels=c(c1, c2, c3))
    df$exprs<-c(as.numeric(t1[g,]), as.numeric(t2[g,]),
                as.numeric(t3[g,]))
    df$std<-c(as.numeric(t1.std[g,]), as.numeric(t2.std[g,]), 
              as.numeric(t3.std[g,]))
    df$ymin<-df$exprs - df$std
    df$ymax<-df$exprs + df$std
    
    p[[r]]<-ggplot(df, aes(x, exprs, group = Condition, colour = Condition)) +
      geom_line(size = 1.5) + theme(axis.text.x = element_blank()) + 
      scale_color_manual(values = col) +
      geom_ribbon(aes(x = x, ymin=ymin, ymax=ymax), linetype="dotted", alpha=0.1) + 
      labs(x = "Pseudotime", y = "Expression") + 
      ggtitle(g) 
    p[[r]]<-p[[r]] + theme(text = element_text(face = "plain", color="black"), 
                           title = element_text(size = 36), 
                           legend.text = element_text(size = 32),
                           axis.text = element_text(size = 36, color="black"),
                           axis.text.x = element_text(angle = 30, hjust = 1, color="black"),
                           axis.line = element_line(colour = "black", size = 1.5),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.35, "cm")) +
            scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())
    r<-r+1
  }
  if(r==2){return(p[[1]])}
  if(r>2){return(plot_grid(plotlist = p))}
}
}
