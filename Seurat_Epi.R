rm(list = ls())
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/seurat流程")

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
exp <- readRDS("Epi_exp.rds")
meta <- read.csv("Epi_meta.csv", row.names = 1)

##初始化Seurat对象
Epi <- CreateSeuratObject(counts = exp, project = "Epi", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = Epi@var.genes, npcs = 20, verbose = FALSE)

##添加tissue信息
tissue <- meta$tissue
names(tissue) <- colnames(x = Epi)
Epi <- AddMetaData(
  object = Epi,
  metadata = tissue,
  col.name = "tissue"
)

p1 <- DimPlot(object = Epi, reduction = "pca", pt.size = .1)
p2 <- VlnPlot(object = Epi, features = "PC_1", pt.size = .1)
plot_grid(p1,p2)

##Run Harmony
Epi <- Epi %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(Epi, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = Epi, reduction = "harmony", pt.size = .1)
p2 <- VlnPlot(object = Epi, features = "harmony_1", pt.size = .1)
plot_grid(p1,p2)

##下游分析
Epi <- Epi %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.6) %>%
  identity()

##可视化
DimPlot(Epi, reduction = "umap", label = TRUE, pt.size = .1)

DimPlot(Epi, reduction = "umap", group.by = "tissue", pt.size = .1)

FeaturePlot(Epi, features = c("CA1","APCDD1","TGFB1","GUCA2B","SMOC2",
                              "CYP2W1","MUC2","REG1A","BST2"),
            pt.size = .1)

##寻找各群的标记
Epi.markers <- FindAllMarkers(Epi,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

write.csv(Epi.markers,"Epi_marker.csv")


##细胞类型注释
new.cluster.ids <- c("Benign","Enterocyte","SPINK4+ CLCA1+ DEFA6+ REG4+",
                     "Benign","Malignant","TUBA1B+ H2AFZ+ HMGB2+ HIST1H4C+",
                     "Malignant","TUBA1B+ H2AFZ+ HMGB2+ HIST1H4C+",
                     "Enterocyte","TUBA1B+ H2AFZ+ HMGB2+ HIST1H4C+",
                     "Goblet","Enterocyte","Malignant","Enterocyte",
                     "Enterocyte","Enteroendocrine")

names(new.cluster.ids) <- levels(Epi)
Epi <- RenameIdents(Epi, new.cluster.ids)

table(Epi@active.ident)
ident <- data.frame(Epi@active.ident)
Epi@meta.data$cell_type <- ident$Epi.active.ident

##绘图
##提取前两个主成分数据
# extact PC ranges
pc12 <- Embeddings(object = Epi,reduction = 'umap') %>%
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
color <- c("#f39c90","#6dba95","#ee787d","#4d62ae","#b37eb7","#53b0d9","#8189b0")

plot1 <- DimPlot(Epi, reduction = 'umap', label = T,
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
cd_genes <- c("GUCA2A","GUCA2B","CA4","CA2",
              "MUC2","TFF3",
              "CHGA","TTR",
              "SPINK4","CLCA1","DEFA6","REG4",
              "REG1A","APCDD1","SELENBP1","SMOC2",
              "TUBA1B","H2AFZ","HMGB2","HIST1H4C",
              "MMP7","LYZ","ERO1A","BST2")

plot2 <- DotPlot(object = Epi, features = cd_genes,
                 cols = c("#FFFFFF", "#BB0021"),
                 col.min = -1,dot.scale = 5)+ 
  theme(axis.text.x = element_text(hjust = 1,angle = 45,size = 12),
        axis.text.y = element_text(size = 12))+
  xlab('')+ylab('')
plot2
ggsave('Epi_dotplot.pdf',plot2,width = 12,height = 5)


##寻找恶性细胞的标记
markers <- FindMarkers(Epi,
                       ident.1 = "Malignant",
                       group.by = "cell_type",
                       min.pct = 0.25,
                       logfc.threshold = 0.6,
                       only.pos = T)

write.csv(markers,"Malignant_DEG.csv")

##保存文件
write.csv(Epi@meta.data, "Epi_seurat_meta.csv")
saveRDS(Epi, file = "Epi_seurat.rds")
