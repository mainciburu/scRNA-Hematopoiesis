library(ggplot2)
library(effectsize)
library(ggsignif)

celltype.y<-read.table("/home/mainciburu/scRNA/scenic/data/young_CellType.txt", stringsAsFactors=F)
celltype.s<-read.table("/home/mainciburu/scRNA/scenic/data/senior_CellType.txt", stringsAsFactors=F)

celltype.levels<-c("HSC", "LMPP", "GMP", "GMP_Granulocytes", 
                    "Monocytes", "pDC", "CLP", "T_NK", "ProB", 
                    "MEP", "Megakaryocytes", "Erythroid_early", 
                    "Erythroid_late", "Basophils")

#######################################
##########   Pseudotime  ##############
#######################################

pst.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/pseudotime.csv", row.names = 1)
pst.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/pseudotime.csv", row.names = 1)

pst.y$Condition<-"Young"
pst.s$Condition<-"Elderly"

df<-rbind(pst.y,pst.s)
colnames(df)[1]<-"Pseudotime"

### Test for distribution
ks<-ks.test(pst.y$X0, pst.s$X0)
a<-ecdf(pst.y$X0)    # cumulative distribution
b<-ecdf(pst.s$X0)

pdf("/home/mainciburu/scRNA/figures/supp_figure4/pseudotime_density.pdf")
ggplot(df, aes(Pseudotime, color = Condition)) + geom_density() + theme_bw() + 
	   scale_color_manual(values = c(Young="#D73027", Elderly ="#4575B4")) + 
	   labs(title = "Pseudotime Distribution",
            subtitle = paste0('KS p-value:', signif(ks$p.val, 3),
                                '; D value:', signif(ks$statistic,3)))
dev.off()

### Test for median
df$CellType<-c(celltype.y$CellType[match(rownames(pst.y), rownames(celltype.y))],
			   celltype.s$CellType[match(rownames(pst.s), rownames(celltype.s))])
df$CellType<-factor(df$CellType, levels = celltype.levels)
df<-df[!is.na(df$CellType),]
df$Condition<-factor(df$Condition, levels = c("Young", "Elderly"))

wx<-c()
rb<-c()
for(ix in unique(df$CellType)){
	print(paste0(ix, "--------"))
	dfi<-df[df$CellType==ix,]
	wxi<-wilcox.test(dfi$Pseudotime[dfi$Condition=="Young"],
		            dfi$Pseudotime[dfi$Condition=="Elderly"])
	wx<-c(wx, wxi$p.val)
	rbi<-rank_biserial(dfi$Pseudotime[dfi$Condition=="Young"],     
		            dfi$Pseudotime[dfi$Condition=="Elderly"])
	rb<-c(rb, rbi$r_rank_biserial)
}
names(wx)<-names(rb)<-unique(df$CellType)

sig<-data.frame(Identity=celltype.levels, 
                 pval=1,
                 significance=factor("NS", levels = c("NS", "*","**", "***")))
i<-match(sig$Identity, names(wx))
sig$pval<-wx[i]
sig$significance[sig$pval<0.05]<-"*"
sig$significance[sig$pval<0.01]<-"**"
sig$significance[sig$pval<0.001]<-"***"

pdf("/home/mainciburu/scRNA/figures/supp_figure4/pseudotime_per_celltype.pdf", width = 10, height = 5)
n<-nrow(sig)
ypos<-1.1
ggplot(df, aes(CellType, Pseudotime, fill = Condition)) + geom_violin(scale = "width") + theme_bw() +
       theme(axis.text.x=element_text(angle = 30, hjust = 1)) + ylim(c(0,1.2)) +
	   scale_fill_manual(values = c(Young="#D73027", Elderly ="#4575B4")) +
	   geom_signif(annotations = sig$significance,
                        y_position = rep(ypos, n),
                        xmin = 1:n-0.3, xmax = 1:n + 0.3, tip_length = 0.01, 
                        size = 1, textsize = 8)
dev.off()


#######################################
####   Differentiation Potential  #####
#######################################

dp.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/diff_potential.csv", row.names = 1)
dp.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/diff_potential.csv", row.names = 1)

ks.test(dp.y$X0, dp.s$X0)
a<-ecdf(dp.y$X0)    # cumulative distribution
b<-ecdf(dp.s$X0)

dp.y$Condition<-"Young"
dp.s$Condition<-"Elderly"

df<-rbind(dp.y,dp.s)
colnames(df)[1]<-"DifferentiationPotential"

pdf("/home/mainciburu/scRNA/figures/supp_figure4/diff_pot_density.pdf")
ggplot(df, aes(DifferentiationPotential, color = Condition)) + geom_density() + theme_bw() + 
	   scale_color_manual(values = c(Young="#D73027", Elderly ="#4575B4")) + 
	   labs(title = "Differentiation Potential Distribution",
            subtitle = paste0('KS p-value:', signif(ks$p.val, 3),
                                '; D value:', signif(ks$statistic,3)))
dev.off()

### Test for median
df$CellType<-c(celltype.y$CellType[match(rownames(dp.y), rownames(celltype.y))],
			   celltype.s$CellType[match(rownames(dp.s), rownames(celltype.s))])
df$CellType<-factor(df$CellType, levels = celltype.levels)
df<-df[!is.na(df$CellType),]
df$Condition<-factor(df$Condition, levels = c("Young", "Elderly"))
wx<-c()
rb<-c()
for(ix in unique(df$CellType)){
	print(paste0(ix, "--------"))
	dfi<-df[df$CellType==ix,]
	wxi<-wilcox.test(dfi$DifferentiationPotential[dfi$Condition=="Young"],
		            dfi$DifferentiationPotential[dfi$Condition=="Elderly"])
	wx<-c(wx, wxi$p.val)
	rbi<-rank_biserial(dfi$DifferentiationPotential[dfi$Condition=="Young"],     
		            dfi$DifferentiationPotential[dfi$Condition=="Elderly"])
	rb<-c(rb, rbi$r_rank_biserial)
}
names(wx)<-names(rb)<-unique(df$CellType)

sig<-data.frame(Identity=celltype.levels, 
                 pval=1,
                 significance=factor("NS", levels = c("NS", "*","**", "***")))
i<-match(sig$Identity, names(wx))
sig$pval<-wx[i]
sig$significance[sig$pval<0.05]<-"*"
sig$significance[sig$pval<0.01]<-"**"
sig$significance[sig$pval<0.001]<-"***"

pdf("/home/mainciburu/scRNA/figures/supp_figure4/diff_pot_per_celltype.pdf", width = 10, height = 5)
n<-nrow(sig)
ypos<-2.1
ggplot(df, aes(CellType, DifferentiationPotential, fill = Condition)) + geom_violin(scale = "width") + theme_bw() +
       theme(axis.text.x=element_text(angle = 30, hjust = 1)) + ylim(c(0,2.3)) +
       labs(y = "Differentiation Potential") +
	   scale_fill_manual(values = c(Young="#D73027", Elderly ="#4575B4")) +
	   geom_signif(annotations = sig$significance,
                        y_position = rep(ypos, n),
                        xmin = 1:n-0.3, xmax = 1:n + 0.3, tip_length = 0.01, 
                        size = 1, textsize = 8)
dev.off()

#######################################
########   Branch Probability  ########
#######################################

bp.y<-read.csv("/home/mainciburu/scRNA/palantir/results/young/branch_probs.csv", row.names = 1)
bp.s<-read.csv("/home/mainciburu/scRNA/palantir/results/senior/branch_probs.csv", row.names = 1)

df.y<-data.frame(Monocytes_BP = bp.y$Monocytes, CellType = celltype.y$CellType[match(rownames(bp.y), rownames(celltype.y))], Condition = "Young")
df.s<-data.frame(Monocytes_BP = bp.s$Monocytes, CellType = celltype.s$CellType[match(rownames(bp.s), rownames(celltype.s))], Condition = "Elderly")
df<-rbind(df.y, df.s)
df$CellType<-factor(df$CellType, levels = celltype.levels)
df<-df[!is.na(df$CellType),]
df$Condition<-factor(df$Condition, levels = c("Young", "Elderly"))
wx<-c()
rb<-c()
for(ix in unique(df$CellType)){
	print(paste0(ix, "--------"))
	dfi<-df[df$CellType==ix,]
	wxi<-wilcox.test(dfi$Monocytes_BP[dfi$Condition=="Young"],
		            dfi$Monocytes_BP[dfi$Condition=="Elderly"])
	wx<-c(wx, wxi$p.val)
	rbi<-rank_biserial(dfi$Monocytes_BP[dfi$Condition=="Young"],     
		            dfi$Monocytes_BP[dfi$Condition=="Elderly"])
	rb<-c(rb, rbi$r_rank_biserial)
}
names(wx)<-names(rb)<-unique(df$CellType)

sig<-data.frame(Identity=celltype.levels, 
                 pval=1,
                 significance=factor("NS", levels = c("NS", "*","**", "***")))
i<-match(sig$Identity, names(wx))
sig$pval<-wx[i]
sig$significance[sig$pval<0.05]<-"*"
sig$significance[sig$pval<0.01]<-"**"
sig$significance[sig$pval<0.001]<-"***"

## Only celltypes in the monocyte branch
df<-df[df$CellType%in%c("HSC", "LMPP", "GMP", "GMP_Granulocytes", "Monocytes"),]
sig<-sig[sig$Identity%in%c("HSC", "LMPP", "GMP", "GMP_Granulocytes", "Monocytes"),]
pdf("/home/mainciburu/scRNA/figures/supp_figure4/monocytes_branch_prob_per_celltype.pdf", width = 8, height = 5)
n<-nrow(sig)
ypos<-1.1
ggplot(df, aes(CellType, Monocytes_BP, fill = Condition)) + geom_violin(scale = "width") + theme_bw() +
       theme(axis.text.x=element_text(angle = 30, hjust = 1)) + ylim(c(0,1.3)) +
       labs(y = "Monocytes Branch Probability") +
	   scale_fill_manual(values = c(Young="#D73027", Elderly ="#4575B4")) +
	   geom_signif(annotations = sig$significance,
                        y_position = rep(ypos, n),
                        xmin = 1:n-0.3, xmax = 1:n + 0.3, tip_length = 0.01, 
                        size = 1, textsize = 8)
dev.off()

