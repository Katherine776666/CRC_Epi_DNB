rm(list = ls())
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle")

##加载R包
{
  library(Matrix)
  library(monocle)
  library(patchwork)
  library(Seurat)
  library(dplyr)
  library(Cairo)
  library(ggplot2)
}

##导入数据
exp <- readRDS("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/seurat流程/Epi_exp.rds")
genemeta <- read.csv("Epi_ac_meta.csv", row.names = 1)
diff_genes <- read.csv("Malignant_DEG.csv", row.names = 1)

index1 <- which(colnames(exp) %in% rownames(genemeta))
genecount <- exp[,index1]
dim(genecount)

##进行monocle对象的构建
RA_matrix <- as(as.matrix(genecount), "sparseMatrix")
class(RA_matrix)

feature_ann <- data.frame(gene_id = rownames(RA_matrix), gene_short_name = rownames(RA_matrix))
rownames(feature_ann) <- rownames(RA_matrix)
RA_fd <- new("AnnotatedDataFrame", data = feature_ann)

sample_ann <- genemeta
rownames(sample_ann) <- colnames(RA_matrix)
RA_pd <- new("AnnotatedDataFrame", data =sample_ann)

##创建monocle对象
RA.cds <- newCellDataSet(RA_matrix,
                         phenoData = RA_pd,
                         featureData = RA_fd,
                         expressionFamily = negbinomial.size())

HSMM <- RA.cds

##归一化
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

##过滤低质量细胞
HSMM <- detectGenes(HSMM, min_expr = 3)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)

##挑选差异显著的基因做降维，给细胞排序
ordering_genes <- row.names(subset(diff_genes, p_val_adj < 0.05))
length(ordering_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)

##降维
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

##保存meta信息和拟时序对象
saveRDS(HSMM,file = "Epi_monocle_Malignant_DEG.rds")
write.csv(pData(HSMM),"Epi_monocleMeta_Malignant_DEG.csv")

##绘图
pData(HSMM)$State <- as.factor(pData(HSMM)$State)
pData(HSMM)$pseudotime <- as.factor(pData(HSMM)$pseudotime)

p1 <- plot_cell_trajectory(HSMM, color_by="Pseudotime")
p1
ggsave("Pseudotime.pdf",p1,width = 8,height = 6)


plot_cell_trajectory(HSMM, color_by="State", show_tree=F, 
                     show_branch_points = F,show_cell_names =F,
                     show_state_number = F,show_backbone = T)

p2 <- plot_cell_trajectory(HSMM, color_by="cell_type", show_tree=F, 
                           show_branch_points = F,show_cell_names =F,
                           show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#f39c90","#4d62ae","#b37eb7"))
p2
ggsave("Pseudotime_celltype.pdf",p2,width = 8,height = 6)


p3 <- plot_cell_trajectory(HSMM, color_by="pseudotime", show_tree=F, 
                           show_branch_points = F,show_cell_names =F,
                           show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#c04327","#dd6f56","#e79a88",
                                 "#c491a5","#8c8ec1","#2a4986"))
p3
ggsave("Pseudotime_after.pdf",p3,width = 8,height = 6)


##制作DNB对象
meta <- read.csv("Epi_monocleMeta_Malignant_DEG.csv",
                 header = T,row.names = 1)

test1 <- meta["pseudotime"]
test1 <- as.data.frame(t(test1))

test1 <- test1[, order(colnames(test1))]
genecount <- genecount[, order(colnames(genecount))]

test2 <- rbind(test1,genecount)
write.csv(test2,"EpiALL_DNBcount_raw(labelled_pseudotime).csv")


##基因表达可视化
genetocheck <- c("REG1A","APCDD1","SELENBP1","MMP7","FABP1","ERO1A",
                 "GAPDH","EEF2","HSPA8","EEF1A1")
for (i in 1:10) {
  gene <- genetocheck[i]
  sub1 <- HSMM[gene,]
  plot_genes_jitter(sub1,grouping = "pseudotime",color_by = "pseudotime",plot_trend = T)+
    scale_color_manual(values = c("#c04327","#dd6f56","#e79a88",
                                  "#c491a5","#8c8ec1","#2a4986"))+
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))+
    theme(legend.position = 'none')+
    xlab('Pseudotime')+ylab('Expression')
  ggsave(paste(gene,'_expression.pdf'),width = 8,height = 6)
}


##功能基因集拟时序热图绘制
###导入基因集文件
geneset <- read.csv("geneset.csv")

genes <- geneset$gene
Time_diff <- differentialGeneTest(HSMM[genes,], cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

p1 <- plot_pseudotime_heatmap(HSMM[Time_genes,],
                              num_clusters = 4,
                              show_rownames=T,
                              return_heatmap=T)

ggsave("heatmap.pdf", p1, width = 8, height = 6)
