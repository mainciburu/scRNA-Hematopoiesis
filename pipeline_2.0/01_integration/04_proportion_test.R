############################
# Test for cell proportions
############################

library(Seurat)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(reshape)
library(xlsx)
source("/home/mainciburu/scRNA/colors.r")


### get cell numbers
young<-readRDS("/home/mainciburu/scRNA/young/seurat_young_v3.rds")
tt.y<-table(young$CellType2, young$Patient)
rm(young)

senior<-readRDS("/home/mainciburu/scRNA/senior/seurat_senior_v4.rds")
tt.s<-table(senior$prediction, senior$Patient)
rm(senior)

mds1<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds1.rds")
tt.m1<-table(mds1$prediction)
rm(mds1)

mds2<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds2.rds")
tt.m2<-table(mds2$prediction)
rm(mds2)

mds3<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds3.rds")
tt.m3<-table(mds3$prediction)
rm(mds3)

mds4<-readRDS("/home/mainciburu/scRNA/MDS_paper/seurat_mds4.rds")
tt.m4<-table(mds4$prediction)
rm(mds4)

tt.y <- rbind(tt.y, "not assigned" = c(0, 0, 0, 0, 0))
tt<-cbind(tt.y, tt.s, MDS1 = tt.m1, MDS3 = tt.m2, MDS5 = tt.m3, MDS10 = tt.m4)
tt<-cbind(tt, Young = rowSums(tt.y), Senior = rowSums(tt.s))
write.table(tt, file = "/home/mainciburu/scRNA/cell_number.txt", quote = F, row.names=T, col.names=T)

### Proportions inside condition
# Test for the 5 young samples
InsideYoung<-c()
for(ii in 1:nrow(tt.y)){
    celltype <- rownames(tt.y)[ii]
    x<-prop.test(x = tt.y[ii,], n = colSums(tt.y))
    print(paste0(celltype, "--------"))
    print(c(Statistic = x$statistic, Pval = x$p.value))
}

# Pairwise comparisons
for(ii in 1:nrow(tt.y)){
    celltype <- rownames(tt.y)[ii]
    x<-pairwise.prop.test(x = tt.y[ii,], 
                 n = colSums(tt.y), 
                 p.adjust.method = "BH")
    print(paste0(celltype, "--------"))
    print(x)
    res<-melt(t(x$p.value))
    res<-res[!is.na(res$value),]   
    if(ii==1){
        WithinYoung<-data.frame(t(res$value))
        colnames(WithinYoung)<-paste(res$X1, res$X2, sep = " vs ")
        rownames(WithinYoung)<-celltype
    }
    else{
        WithinYoung<-rbind(WithinYoung, res$value)
        rownames(WithinYoung)[[ii]]<-celltype
    }
}


# Test for the 3 senior samples
for(ii in 1:nrow(tt.s)){
    celltype <- rownames(tt.s)[ii]
    x<-prop.test(x = tt.s[ii,], n = colSums(tt.s))
    print(paste0(celltype, "--------"))
    print(c(Statistic = x$statistic, Pval = x$p.value))
}

# Pairwise comparisons
for(ii in 1:nrow(tt.s)){
    celltype <- rownames(tt.s)[ii]
    x<-pairwise.prop.test(x = tt.s[ii,], 
                 n = colSums(tt.s), 
                 p.adjust.method = "BH")
    print(paste0(celltype, "--------"))
    print(x)
    res<-melt(t(x$p.value))
    res<-res[!is.na(res$value),]   
    if(ii==1){
        WithinSenior<-data.frame(t(res$value))
        colnames(WithinSenior)<-paste(res$X1, res$X2, sep = " vs ")
        rownames(WithinSenior)<-celltype
    }
    else{
        WithinSenior<-rbind(WithinSenior, res$value)
        rownames(WithinSenior)[[ii]]<-celltype
    }
}


### Proportions between young, elderly and MDS patients
# Test for any difference
for(ii in 1:nrow(tt)){
    celltype <- rownames(tt)[ii]
    x<-prop.test(x = tt[ii,c("Young", "Senior", "MDS1", "MDS2", "MDS3", "MDS4")], 
                 n = colSums(tt[,c("Young", "Senior", "MDS1", "MDS2", "MDS3", "MDS4")]))
    print(paste0(celltype, "--------"))
    print(c(Statistic = x$statistic, Pval = x$p.value))
}

# Pairwise comparisons
for(ii in 1:nrow(tt)){
    celltype <- rownames(tt)[ii]
    x<-pairwise.prop.test(x = tt[ii,c("Young", "Senior", "MDS1", "MDS2", "MDS3", "MDS4")], 
                 n = colSums(tt[,c("Young", "Senior", "MDS1", "MDS2", "MDS3", "MDS4")]), 
                 p.adjust.method = "BH")
    print(paste0(celltype, "--------"))
    print(x)
    res<-melt(t(x$p.value))
    res<-res[!is.na(res$value),]   
    if(ii==1){
        BetweenCondition<-data.frame(t(res$value))
        colnames(BetweenCondition)<-paste(res$X1, res$X2, sep = " vs ")
        rownames(BetweenCondition)<-celltype
    }
    else{
        BetweenCondition<-rbind(BetweenCondition, res$value)
        rownames(BetweenCondition)[[ii]]<-celltype
    }
}


### Write results
# remove not assigned row
WithinSenior<-WithinSenior[-15,]
BetweenCondition<-BetweenCondition[-15,]
ProportionTable<-cbind(WithinYoung, WithinSenior, BetweenCondition)
write.xlsx(ProportionTable, file = "scRNA/proportion_test_results.xlsx")



