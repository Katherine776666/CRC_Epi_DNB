rm(list = ls())
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277")

##加载R包
{
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(Cairo)
  library(ggplot2)
  library(cowplot)
  library(harmony)
}

##导入数据
rawcount <- readRDS("GSE161277_exp.rds")

##初始化Seurat对象
CRC <- CreateSeuratObject(counts = rawcount, project = "GSE161277", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = CRC@var.genes, npcs = 20, verbose = FALSE)

##线性降维后绘图
p1 <- DimPlot(object = CRC, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = CRC, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)

##Run Harmony
CRC <- CRC %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(CRC, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = CRC, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = CRC, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)

##下游分析
CRC <- CRC %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  identity()

##可视化
DimPlot(CRC, reduction = "umap", label = TRUE, pt.size = .1)

FeaturePlot(CRC,features = c("GUCA2A","MUC2","CHGA","EPCAM","CD3D","CD8A",
                             "KLRD1","MS4A1","MZB1","CD68","CD163","DCN"))


##寻找各群的标记
CRC.markers <- FindAllMarkers(CRC,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

write.csv(CRC.markers, "CRC_markers.csv")

##细胞类型注释
new.cluster.ids <- c("T cell","T cell","Epithelial cell","Follicular B cell",
                     "Epithelial cell","Epithelial cell","T cell",
                     "Epithelial cell","Follicular B cell","Plasma B cell",
                     "Epithelial cell","Macrophage","Macrophage",
                     "Epithelial cell","Macrophage","Epithelial cell",
                     "Epithelial cell","Epithelial cell","Epithelial cell",
                     "Follicular B cell","Epithelial cell","Epithelial cell",
                     "Fibroblast","Epithelial cell")


names(new.cluster.ids) <- levels(CRC)
CRC <- RenameIdents(CRC, new.cluster.ids)

ident <- data.frame(CRC@active.ident)
CRC@meta.data$cell_type <- ident$CRC.active.ident
table(CRC@meta.data$cell_type)

##绘图
##提取前两个主成分数据
# extact PC ranges
pc12 <- Embeddings(object = CRC,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)

##构造坐标轴需要的标签和位置信息
# get botomn-left coord
lower <- floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2))) - 2

# get relative line length
linelen <- abs(0.3*lower) + lower

# mid point
mid <- abs(0.3*lower)/2 + lower

# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('UMAP_2','UMAP_1'),each = 2))

# axies label
label <- data.frame(lab = c('UMAP_2','UMAP_1'),angle = c(90,0),
                    x = c(lower - 3,mid),y = c(mid,lower - 2.5))

##可视化
# plot
color <- c("#6dba95","#ee787d","#b37eb7","#53b0d9","#4d62ae","#8189b0")

plot1 <- DimPlot(CRC, reduction = 'umap', label = T,
                 pt.size = .1, cols = color) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot1
ggsave('Celltype_umap.pdf',plot1,width = 8,height = 8)

##绘制Dotplot
cd_genes <- c("CD3D","CD3E","CD3G","KLRB1","KLRD1",
              "EPCAM","KRT18","KRT8","TFF3","KRT19",
              "CD79A","CD79B","MS4A1","LTB","MZB1",
              "CD68","CD14","CD86","CD163","IL1B",
              "COL1A1","DCN","COL1A2","COL3A1","COL6A2")

plot2 <- DotPlot(object = CRC, features = cd_genes,
                 cols = c("#FFFFFF", "#BB0021"),
                 col.min = -1,dot.scale = 5)+ 
  theme(axis.text.x = element_text(hjust = 1,angle = 45,size = 12),
        axis.text.y = element_text(size = 12))+
  xlab('')+ylab('')
plot2
ggsave('CRC_dotplot.pdf',plot2,width = 12,height = 5)


##保存文件
write.csv(CRC@meta.data, "meta.csv")
saveRDS(CRC, file = "GSE161277_seurat.rds")

Epi <- subset(x=CRC, cell_type == "Epithelial cell")
Epi_exp <- as.matrix(Epi@assays$RNA@counts)
saveRDS(Epi_exp, file = "Epi_exp.rds")

