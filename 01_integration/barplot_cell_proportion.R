### Barplot
# **********
# Here I am using the proportion mean and sd resulting from the 5 young and 3 seniors
# For the proportion test, I didn't use this. I directly use cell numbers from the integrated data set
# **********

# Function to calculate the mean and the standard deviation
# for each group
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

prop.tt<-prop.table(tt, margin = 2)[,c(11,12)]   # don't include Young and Senior // not assigned

df<-melt(t(prop.tt))
colnames(df)<-c("Patient", "CellType", "Proportion")
df$Condition<-as.character(df$Patient)
df$Condition[grep("young", df$Patient)]<-"Young"
df$Condition[grep("senior", df$Patient)]<-"Elderly"
df$Condition<-factor(df$Condition, levels = c("Young", "Elderly", "MDS", "AML"))
df2<-data_summary(df, varname = "Proportion", groupnames = c("Condition", "CellType"))
df2$sd[is.na(df2$sd)]<-0
df2$CellType<-factor(df2$CellType, levels = rownames(prop.tt))


pdf("/home/mainciburu/scRNA/figures/cell_proportion_grouped_barplot.pdf", width = 10, useDingbats = F)
ggplot(df2, aes(CellType, Proportion, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = Proportion - sd, ymax = Proportion + sd), width=.2, position=position_dodge(.9)) + 
  theme_minimal() +
  scale_fill_manual(values = col.condition) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(text = element_text(face = "bold", color = "black"), 
        axis.text = element_text(size = 18, color = "black"), 
        axis.text.x = element_text(angle = 30, hjust = 1), 
        axis.title = element_text(size = 20), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) + theme(legend.position = c(0.8, 0.8))
dev.off()
