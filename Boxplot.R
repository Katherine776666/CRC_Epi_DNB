rm(list = ls()) 
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/基因调控路线")

##加载R包
{
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  library(knitr)
  library(forcats)
}

##导入数据
Epi <- readRDS("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/seurat流程/Epi_seurat.rds")
exp <- as.matrix(Epi@assays$RNA@data)

group <- read.csv("group.csv", row.names = 1)
exp <- data.frame(t(exp))

index1 <- which(rownames(exp) %in% rownames(group))
exp <- exp[index1,]

expr1 <- exp %>%
  select(., FOS,JUN,TP53)
expr1$group <- group$group

##FOS表达箱线图
plot1 <- ggplot(aes(x = group, y = FOS,fill = group),data = expr1) +
  stat_boxplot(geom = 'errorbar',width=0.25,cex=0.5)+
  geom_boxplot(outlier.size = -1,width=0.5,fill=c("#6c6cb5","#e07979")) +
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name = "") +
  ggtitle("Expression of FOS") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        text = element_text(size = 10),
        axis.text.x=element_text(size = 10,angle = 0),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  geom_signif(comparisons = list(c("Before-deterioration stage","Deteriorated stage")),
              test="wilcox.test",
              tip_length = 0.03,size = 0.8,textsize = 4)
plot1
ggsave("Expression of FOS.pdf",plot1,width = 3,height = 6)

##JUN表达箱线图
plot2 <- ggplot(aes(x = group, y = JUN,fill = group),data = expr1) +
  stat_boxplot(geom = 'errorbar',width=0.25,cex=0.5)+
  geom_boxplot(outlier.size = -1,width=0.5,fill=c("#6c6cb5","#e07979")) +
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name = "") +
  ggtitle("Expression of JUN") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        text = element_text(size = 10),
        axis.text.x=element_text(size = 10,angle = 0),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  geom_signif(comparisons = list(c("Before-deterioration stage","Deteriorated stage")),
              test="wilcox.test",
              tip_length = 0.03,size = 0.8,textsize = 4)
plot2
ggsave("Expression of JUN.pdf",plot2,width = 3,height = 6)

##TP53表达箱线图
plot3 <- ggplot(aes(x = group, y = TP53,fill = group),data = expr1) +
  stat_boxplot(geom = 'errorbar',width=0.25,cex=0.5)+
  geom_boxplot(outlier.size = -1,width=0.5,fill=c("#6c6cb5","#e07979")) +
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name = "") +
  ggtitle("Expression of TP53") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        text = element_text(size = 10),
        axis.text.x=element_text(size = 10,angle = 0),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  geom_signif(comparisons = list(c("Before-deterioration stage","Deteriorated stage")),
              test="wilcox.test",
              tip_length = 0.03,size = 0.8,textsize = 4)
plot3
ggsave("Expression of TP53.pdf",plot3,width = 3,height = 6)


##凋亡相关基因表达箱线图
expr2 <- exp %>% select(., BBC3,CASP6,GADD45A,XIAP)
expr2$group <- group$group
expr2 <- tibble::rownames_to_column(expr2, var = "Cell_id")
mydata <- melt(expr2,
               id.vars=c("Cell_id", "group"),
               variable.name="gene",value.name="expression")
mydata$gene <- fct_inorder(mydata$gene)

compaired <- list(c("Before-deterioration stage","Deteriorated stage"))

plot4 <- ggplot(aes(x = gene, y = expression,fill = group),data = mydata) +
  stat_boxplot(geom = 'errorbar',width=0.25,cex=0.5)+
  geom_boxplot(outlier.size = -1,width=0.5) +
  scale_fill_manual(values = c("#6c6cb5","#e07979"))+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name = "") +
  ggtitle("Expression of apoptosis-related genes") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        text = element_text(size = 10),
        axis.text.x=element_text(size = 10,angle = 0),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  geom_signif(comparisons = compaired,test="wilcox.test",
              map_signif_level = T,step_increase = 0,na.rm = FALSE,
              tip_length = 0.03,size = 0.8,textsize = 4)
plot4
ggsave("Expression of apoptosis-related genes.pdf",plot4,width = 10,height = 6)

