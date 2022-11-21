rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/GSVA")

##加载R包
{
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
}

##导入数据
exp <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle/EpiALL_DNBcount_raw(labelled_pseudotime).csv",
                header = T, row.names = 1, check.names = F)
geneset <- read.csv("geneset.csv", header = T)

##构建GMT文件
name <- unique(geneset$term)
description <- rep(NA,length(name))
names(description) <- name

genes <- lapply(name, function(name){
  as.vector(geneset[geneset$term == name,"gene"])
})
names(genes) <- name

gmtinput <- list(name=name,description=description,genes=genes)

get_gmt <- function(gmtinput,filename){
  output <- file(filename, open="wt")
  lapply(gmtinput[["name"]],function(name){
    outlines = paste0(c(name, gmtinput[["description"]][[name]],
                        gmtinput[["genes"]][[name]]),collapse='\t')
    writeLines(outlines, con=output)
  })
  close(output)
}
get_gmt(gmtinput=gmtinput,filename="geneset.gmt")


##导入GSVA分析所需文件
geneSets <- getGmt("geneset.gmt")

mydata <- as.matrix(exp)
mydata <- round(mydata)

##GSVA富集分析(这一步比较耗时)
Es <- gsva(mydata, geneSets, parallel.sz=1,kcdf="Poisson")
write.csv(Es, "Es.csv")


##limma进行差异通路分析
Es <- read.csv("Es.csv", row.names = 1, check.names = F)

###导入分组信息(A-C)
group <- read.csv("group_AVSC.csv", row.names = 1)

index1 <- which(colnames(Es) %in% rownames(group))
Es_1 <- Es[,index1]

###分组矩阵
design.matrix <- model.matrix(~0+factor(group$group))
colnames(design.matrix) <- c('After','Critical')
rownames(design.matrix) <- rownames(group)

###比较矩阵
contracts.matrix <- makeContrasts(After-Critical,
                                  levels=design.matrix)

###差异分析
fit1 <- lmFit(Es_1, design.matrix)
fit2 <- contrasts.fit(fit1, contracts.matrix)
fit3 <- eBayes(fit2)
summary(decideTests(fit3,adj.P.Val = 0.05))

###结果输出
pathway_limma <- topTable(fit3,
                          number = Inf,
                          sort.by = 'logFC') %>% filter(.,adj.P.Val<0.05)
pathway_limma <- pathway_limma[order(pathway_limma$logFC,decreasing = T),]
write.csv(pathway_limma, "pathway_limma_AC.csv")


###导入分组信息(B-C)
group2 <- read.csv("group_BVSC.csv", row.names = 1)

index2 <- which(colnames(Es) %in% rownames(group2))
Es_2 <- Es[,index2]

###分组矩阵
design.matrix_2 <- model.matrix(~0+factor(group2$group))
colnames(design.matrix_2) <- c('Before','Critical')
rownames(design.matrix_2) <- rownames(group2)

###比较矩阵
contracts.matrix_2 <- makeContrasts(Before-Critical,
                                    levels=design.matrix_2)

###差异分析
fit4 <- lmFit(Es_2, design.matrix_2)
fit5 <- contrasts.fit(fit4, contracts.matrix_2)
fit6 <- eBayes(fit5)
summary(decideTests(fit6,adj.P.Val = 0.05))

###结果输出
pathway_limma_2 <- topTable(fit6,
                            number = Inf,
                            sort.by = 'logFC') %>% filter(.,adj.P.Val<0.05)
pathway_limma_2 <- pathway_limma_2[order(pathway_limma_2$logFC,decreasing = T),]
write.csv(pathway_limma_2, "pathway_limma_BC.csv")


##绘图数据整理
Es <- read.csv("Es.csv", row.names = 1, check.names = F)
colnames(Es) <- exp[1,]

##分群均值处理
index1 <- names(Es) ==1
index2 <- names(Es) ==2
index3 <- names(Es) ==3
index4 <- names(Es) ==4
index5 <- names(Es) ==5
index6 <- names(Es) ==6

exp_1 <- Es[,index1]
exp_2 <- Es[,index2]
exp_3 <- Es[,index3]
exp_4 <- Es[,index4]
exp_5 <- Es[,index5]
exp_6 <- Es[,index6]

Cluster1 <- apply(exp_1,1,mean)
Cluster2 <- apply(exp_2,1,mean)
Cluster3 <- apply(exp_3,1,mean)
Cluster4 <- apply(exp_4,1,mean)
Cluster5 <- apply(exp_5,1,mean)
Cluster6 <- apply(exp_6,1,mean)

mean_Es <- data.frame(Cluster1,Cluster2,Cluster3,
                      Cluster4,Cluster5,Cluster6)

##绘图
##构建列注释信息
annotation_col <- data.frame(
  Celltype = factor(rep(c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell",
                          "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell",
                          "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell",
                          "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell",
                          "LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell",
                          "MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell"),1)))

rownames(annotation_col) <- colnames(mean_Es)

###自定义注释信息的颜色列表
anno_colors <- list(
  Celltype = c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell" = "#c04327",
               "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell" = "#dd6f56",
               "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell" = "#e79a88",
               "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell" = "#c491a5",
               "LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell" = "#8c8ec1",
               "MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell" = "#2a4986"))

###数据scale标准化处理
Es_scale <- mean_Es %>% t() %>% scale() %>% t()

plot <- pheatmap(Es_scale,cluster_cols = F,
                 cluster_rows = T, border_color = "white",
                 color = colorRampPalette(c("#3758a9","white","#BC2A33"))(100),
                 annotation_col = annotation_col,
                 annotation_colors = anno_colors,
                 show_colnames =F, show_rownames = T,
                 treeheight_row = 10,
                 fontsize_number = 12,
                 fontsize_row = 10,
                 fontsize_col = 10,
                 cellwidth = 25,
                 cellheight = 15)

ggsave("ES_heatmap.pdf",plot,width = 10,height = 6)
