rm(list = ls())
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/FindMarkers")

##加载R包
{
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(Cairo)
  library(ggplot2)
  library(cowplot)
  library(harmony)
  library(dplyr)
  library(pheatmap)
}

##导入数据
exp <- readRDS("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/seurat流程/Epi_exp.rds")
meta <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle/Epi_ac_meta.csv",
                 header = T,row.names = 1)

index1 <- which(colnames(exp) %in% rownames(meta))
exp <- exp[,index1]
dim(exp)

##初始化Seurat对象
EpiAC <- CreateSeuratObject(counts = exp, project = "EpiAC") %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = EpiAC@var.genes, npcs = 20, verbose = FALSE)

#添加pseudotime信息
pseudotime <- meta$pseudotime
names(pseudotime) <- colnames(x = EpiAC)
EpiAC <- AddMetaData(
  object = EpiAC,
  metadata = pseudotime,
  col.name = "pseudotime"
)

##线性降维后绘图
p1 <- DimPlot(object = EpiAC, reduction = "pca", pt.size = .1)
p2 <- VlnPlot(object = EpiAC, features = "PC_1", pt.size = .1)
plot_grid(p1,p2)

##Run Harmony
EpiAC <- EpiAC %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(EpiAC, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = EpiAC, reduction = "harmony", pt.size = .1)
p2 <- VlnPlot(object = EpiAC, features = "harmony_1", pt.size = .1)
plot_grid(p1,p2)

##下游分析
EpiAC <- EpiAC %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.6) %>%
  identity()

##可视化
DimPlot(EpiAC, reduction = "umap", label = T, pt.size = .1)

##寻找拟时序各个时间点的差异表达基因
deg1 <- FindMarkers(EpiAC,
                    ident.1 = 1,
                    group.by = "pseudotime",
                    min.pct = 0.25,
                    logfc.threshold = 0.25)

deg2 <- FindMarkers(EpiAC,
                    ident.1 = 2,
                    group.by = "pseudotime",
                    min.pct = 0.25,
                    logfc.threshold = 0.25)

deg3 <- FindMarkers(EpiAC,
                    ident.1 = 3,
                    group.by = "pseudotime",
                    min.pct = 0.25,
                    logfc.threshold = 0.25)

deg4 <- FindMarkers(EpiAC,
                    ident.1 = 4,
                    group.by = "pseudotime",
                    min.pct = 0.25,
                    logfc.threshold = 0.25)

deg5 <- FindMarkers(EpiAC,
                    ident.1 = 5,
                    group.by = "pseudotime",
                    min.pct = 0.25,
                    logfc.threshold = 0.25)

deg6 <- FindMarkers(EpiAC,
                    ident.1 = 6,
                    group.by = "pseudotime",
                    min.pct = 0.25,
                    logfc.threshold = 0.25)

##保存文件
write.csv(deg1, "cluster1.markers.csv")
write.csv(deg2, "cluster2.markers.csv")
write.csv(deg3, "cluster3.markers.csv")
write.csv(deg4, "cluster4.markers.csv")
write.csv(deg5, "cluster5.markers.csv")
write.csv(deg6, "cluster6.markers.csv")


##绘图数据处理
deg_all <- rbind(deg1,deg2,deg3,deg4,deg5,deg6)
expr <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle/EpiALL_DNBcount_raw(labelled_pseudotime).csv",
                header = T, row.names = 1)
names(expr) <- expr[1,]
expr <- expr[-1,]

index_deg <- which(rownames(expr) %in% rownames(deg_all))
expr_deg <- expr[index_deg,]

##分群均值处理
index1 <- names(expr_deg) ==1
index2 <- names(expr_deg) ==2
index3 <- names(expr_deg) ==3
index4 <- names(expr_deg) ==4
index5 <- names(expr_deg) ==5
index6 <- names(expr_deg) ==6

exp_1 <- expr_deg[,index1]
exp_2 <- expr_deg[,index2]
exp_3 <- expr_deg[,index3]
exp_4 <- expr_deg[,index4]
exp_5 <- expr_deg[,index5]
exp_6 <- expr_deg[,index6]

Cluster1 <- apply(exp_1,1,mean)
Cluster2 <- apply(exp_2,1,mean)
Cluster3 <- apply(exp_3,1,mean)
Cluster4 <- apply(exp_4,1,mean)
Cluster5 <- apply(exp_5,1,mean)
Cluster6 <- apply(exp_6,1,mean)

deg_mathix <- data.frame(Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6)

##绘制热图
###构建列注释信息
annotation_col <- data.frame(
  Celltype = factor(rep(c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell",
                          "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell",
                          "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell",
                          "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell",
                          "LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell",
                          "MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell"),1)))
rownames(annotation_col) <- colnames(deg_mathix)

###自定义注释信息的颜色列表
anno_colors <- list(
  Celltype = c("REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell" = "#c04327",
               "SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell" = "#dd6f56",
               "PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell" = "#e79a88",
               "FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell" = "#c491a5",
               "LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell" = "#8c8ec1",
               "MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell" = "#2a4986"))

###数据scale标准化处理
test <- deg_mathix[apply(deg_mathix, 1, function(x) sd(x)!=0),]
deg_mathix_scale <- test %>% t() %>% scale() %>% t()

plot <- pheatmap(deg_mathix_scale, cluster_cols = F,
                 cluster_rows = T, border_color = "white",
                 color = colorRampPalette(c("#3758a9","white","#BC2A33"))(100),
                 annotation_col = annotation_col,
                 annotation_colors = anno_colors,
                 show_colnames =F, show_rownames = F,
                 treeheight_row = 10,
                 fontsize_number = 18,
                 fontsize_row = 10,
                 fontsize_col = 10,
                 cellwidth = 35)
ggsave("DEG_heatmap.pdf",plot,width = 10,height = 6)
