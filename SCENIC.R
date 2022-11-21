rm(list = ls())
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/SCENIC")

##加载R包
{
  library(Seurat)
  library(dplyr)
  library(SCENIC)
  library(stringr)
  library(limma)
  library(pheatmap)
  library(ggplot2)
}

##导入数据
meta <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle/Epi_ac_meta.csv",
                 header = T,row.names = 1)
Epi <- readRDS("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/seurat流程/Epi_seurat.rds")

##获取表达矩阵
exprMat <- as.matrix(Epi@assays$RNA@data)
index1 <- which(colnames(exprMat) %in% rownames(meta))
exprMat <- exprMat[,index1]

##获取临床表型
cellInfo <- meta[,c(8,2,3)]
colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')

##SCENIC分析
###初始化设置
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=1) 
saveRDS(scenicOptions, file="int/scenicOptions.rds") 

###将表达矩阵中未在cisTarget database中收录的基因去除
genesKept <- geneFiltering(exprMat, scenicOptions)
datExpr_filtered <- exprMat[genesKept, ]
dim(datExpr_filtered)

###共表达网络构建（GENIE3）
##相关性计算
runCorrelation(datExpr_filtered,scenicOptions)

##鉴定潜在靶点
runGenie3(datExpr_filtered, scenicOptions)

##推断共表达模块
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

##推断转录调控网络（regulon）
coexMethod <- c("w001","w005","top50","top5perTarget","top10perTarget","top50perTarget")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethods = coexMethod)

##regulon活性评分与可视化
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions,exprMat = datExpr_filtered)

##AUC结果二进制转换
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions,exprMat = datExpr_filtered)

###保存文件
saveRDS(scenicOptions, file="int/scenicOptions.rds") 

##regulon AUC去除extended TF
regulonAUC <- readRDS("./int/3.4_regulonAUC.Rds")
AUCmatrix <- regulonAUC@assays@data@listData$AUC
AUCmatrix_filter <- AUCmatrix[!str_detect(rownames(AUCmatrix),'[extended]'),]


##limma差异分析(A_C)
###导入分组信息
group1 <- read.csv("group_AVSC.csv", row.names = 1)

index1 <- which(colnames(AUCmatrix_filter) %in% rownames(group1))
AUCmatrix_filter_1 <- AUCmatrix_filter[,index1]

###分组矩阵
design.matrix_1 <- model.matrix(~0+factor(group1$group))
colnames(design.matrix_1) <- c('After','Critical')
rownames(design.matrix_1) <- rownames(group1)

###比较矩阵
contracts.matrix_1 <- makeContrasts(After-Critical,
                                    levels=design.matrix_1)

###差异分析
fit1 <- lmFit(AUCmatrix_filter_1, design.matrix_1)
fit2 <- contrasts.fit(fit1, contracts.matrix_1)
fit3 <- eBayes(fit2)
summary(decideTests(fit3,adj.P.Val = 0.05))

###结果输出
auc_limma_1 <- topTable(fit3,
                        number = Inf,
                        sort.by = 'logFC') %>% filter(.,adj.P.Val<0.05)
auc_limma_1 <- auc_limma_1[order(auc_limma_1$logFC,decreasing = T),]
write.csv(auc_limma_1, "auc_limma_AC.csv")


##limma差异分析(B_C)
###导入分组信息
group2 <- read.csv("group_BVSC.csv", row.names = 1)

index2 <- which(colnames(AUCmatrix_filter) %in% rownames(group2))
AUCmatrix_filter_2 <- AUCmatrix_filter[,index2]

###分组矩阵
design.matrix_2 <- model.matrix(~0+factor(group2$group))
colnames(design.matrix_2) <- c('Before','Critical')
rownames(design.matrix_2) <- rownames(group2)

###比较矩阵
contracts.matrix_2 <- makeContrasts(Before-Critical,
                                    levels=design.matrix_2)

###差异分析
fit4 <- lmFit(AUCmatrix_filter_2, design.matrix_2)
fit5 <- contrasts.fit(fit4, contracts.matrix_2)
fit6 <- eBayes(fit5)
summary(decideTests(fit6,adj.P.Val = 0.05))

###结果输出
auc_limma_2 <- topTable(fit6,
                        number = Inf,
                        sort.by = 'logFC') %>% filter(.,adj.P.Val<0.05)
auc_limma_2 <- auc_limma_2[order(auc_limma_2$logFC,decreasing = T),]
write.csv(auc_limma_2, "auc_limma_BC.csv")


##绘制热图
exp <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle/EpiALL_DNBcount_raw(labelled_pseudotime).csv",
                header = T, row.names = 1,check.names = F)
TF_select <- read.csv("TF_select.csv")
index_TF <- which(rownames(AUCmatrix_filter) %in% TF_select$TF)
AUCmatrix_select <- AUCmatrix_filter[index_TF,]
colnames(AUCmatrix_select) <- exp[1,]

###分群均值处理
index1 <- colnames(AUCmatrix_select) ==1
index2 <- colnames(AUCmatrix_select) ==2
index3 <- colnames(AUCmatrix_select) ==3
index4 <- colnames(AUCmatrix_select) ==4
index5 <- colnames(AUCmatrix_select) ==5
index6 <- colnames(AUCmatrix_select) ==6

AUCmatrix_1 <- AUCmatrix_select[,index1]
AUCmatrix_2 <- AUCmatrix_select[,index2]
AUCmatrix_3 <- AUCmatrix_select[,index3]
AUCmatrix_4 <- AUCmatrix_select[,index4]
AUCmatrix_5 <- AUCmatrix_select[,index5]
AUCmatrix_6 <- AUCmatrix_select[,index6]

Cluster1 <- apply(AUCmatrix_1,1,mean)
Cluster2 <- apply(AUCmatrix_2,1,mean)
Cluster3 <- apply(AUCmatrix_3,1,mean)
Cluster4 <- apply(AUCmatrix_4,1,mean)
Cluster5 <- apply(AUCmatrix_5,1,mean)
Cluster6 <- apply(AUCmatrix_6,1,mean)

mathix <- data.frame(Cluster1,Cluster2,Cluster3,
                     Cluster4,Cluster5,Cluster6)

colnames(mathix) <- c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell",
                      "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell",
                      "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell",
                      "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell",
                      "LYZ+ CST3+ STMN1+ CYP2W1+ Malignant Cell",
                      "MMP7+ FABP1+ TFF1+ CKB+ Malignant Cell")

###构建列注释信息
annotation_col <- data.frame(
  Celltype = factor(rep(c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell",
                          "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell",
                          "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell",
                          "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell",
                          "LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell",
                          "MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell"),1)))

rownames(annotation_col) <- colnames(mathix)

###自定义注释信息的颜色列表
anno_colors <- list(
  Celltype = c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell" = "#c04327",
               "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell" = "#dd6f56",
               "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell" = "#e79a88",
               "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell" = "#c491a5",
               "LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell" = "#8c8ec1",
               "MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell" = "#2a4986"))

###数据scale标准化处理
mathix_scale <- mathix %>% t() %>% scale() %>% t()

plot <- pheatmap(mathix_scale,cluster_cols = F,
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

ggsave("TF_heatmap.pdf",plot,width = 10,height = 6)
