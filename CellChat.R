rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/CellChat")

##加载R包
{
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
}

# Part1 数据输入和处理以及CellChat对象的初始化
## 导入数据
CRC <- readRDS("GSE161277_seurat.rds")
meta <- read.csv("meta_ALL.csv", header = T, row.names = 1)

exp <- CRC@assays$RNA@data

index1 <- which(colnames(exp) %in% rownames(meta))
data.input <- exp[,index1]
dim(data.input)

## 创建CellChat对象
cellchat <- createCellChat(object = data.input,
                           meta = meta,
                           group.by = "labels")

## 导入配体受体相互作用数据库
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)# Show the structure of the database
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB

## 数据预处理
cellchat <- subsetData(cellchat)# subset the expression data of signaling genes for saving computation cost
#future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)# project gene expression data onto PPI network (optional)

# Part2 细胞间通讯网络推断
## 计算通讯概率并推断网络
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

## 提取推断得到的细胞通讯网络
df.net <- subsetCommunication(cellchat)
#df.net <- subsetCommunication(all_cellchat, signaling = 'WNT')#netp返回通路级别数据/sources.use返回信号起点终点/signaling 返回特定信号通路

## 在信号通路水平上推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)

## 计算整合网络
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

# mat <- cellchat@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

# 保存CellChat对象
saveRDS(cellchat, file = "cellchat_All.rds")

# 可视化
## 查看重要通信的信号通路
cellchat@netP$pathways

color <- c("#c04327","#dd6f56","#e79a88","#c491a5","#8c8ec1","#2a4986",
           "#8189b0","#b37eb7","#4d62ae","#53b0d9","#6dba95")
target <- c('4_FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell')
source <- c('1_REG1A+ LEFTY1+ SMOC2+ PCCA+ Epithelial Cell',
            '2_SELENBP1+ APCDD1+ DPEP1+ ID3+ Epithelial Cell',
            '3_PIGR+ CD74+ SPINK4+ LCN2+ Epithelial Cell',
            '4_FABP5+ S100P+ PLA2G2A+ TUBA1B+ Epithelial Cell',
            '5_LYZ+ CST3+ STMN1+ CYP2W1+ Epithelial Cell',
            '6_MMP7+ FABP1+ TFF1+ CKB+ Epithelial Cell',
            'Fibroblast','Follicular B cell',
            'Macrophage','Plasma B cell','T cell')

## 挑选用于可视化的信号通路
GDF <- c("GDF")
GRN <- c("GRN")
CEACAM <- c("CEACAM")

## 弦图
plot1 <- netVisual_chord_gene(cellchat, big.gap = 20,small.gap = 5,
                              targets.use = target,transparency = 0.5,
                              signaling = GDF,show.legend = T,color.use = color)
plot2 <- netVisual_chord_gene(cellchat, big.gap = 20,small.gap = 5,
                              targets.use = target,transparency = 0.5,
                              signaling = GRN,show.legend = T,color.use = color)
plot3 <- netVisual_chord_gene(cellchat, big.gap = 20,small.gap = 5,
                              targets.use = target,transparency = 0.5,
                              signaling = CEACAM,show.legend = T,color.use = color)

pdf('chord.pdf', height = 10,width = 8)
plot1
plot2
plot3
dev.off()

## 热图
plot4 <- netVisual_heatmap(cellchat, color.heatmap = 'Reds',color.use = color,
                           signaling = GDF,targets.use = target,
                           sources.use = source[-c(7,8,9,10,11)],
                           title.name = 'GDF signaling pathway interaction')
plot5 <- netVisual_heatmap(cellchat, color.heatmap = 'Reds',color.use = color,
                           signaling = GRN,targets.use = target,
                           title.name = 'GRN signaling pathway interaction')
plot6 <- netVisual_heatmap(cellchat, color.heatmap = 'Reds',color.use = color,
                           signaling = CEACAM,targets.use = target,
                           title.name = 'CEACAM signaling pathway interaction')
pdf('heatmap.pdf', height = 8,width = 8)
plot4
plot5
plot6
dev.off()

## 点图
plot7 <- netVisual_bubble(cellchat,
                          signaling = "GDF",targets.use =  target, 
                          angle.x = 45, remove.isolate =T,
                          title.name = 'GDF signaling pathway')
plot8 <- netVisual_bubble(cellchat,
                          signaling = "GRN",targets.use =  target, 
                          angle.x = 45, remove.isolate =T,
                          title.name = 'GRN signaling pathway')
plot9 <- netVisual_bubble(cellchat,
                          signaling = "CEACAM",targets.use =  target, 
                          angle.x = 45, remove.isolate =T,
                          title.name = 'CEACAM signaling pathway')
pdf('dot.pdf', height = 8,width = 8)
plot7
plot8
plot9
dev.off()

## 导出prob数据用于标注在热图上
GDF_prob <- cellchat@netP[["prob"]][,,'GDF'] 
GRN_prob <- cellchat@netP[["prob"]][,,'GRN'] 
CEACAM_prob <- cellchat@netP[["prob"]][,,'CEACAM'] 

list_of_datasets <- list('GDF' = GDF_prob,
                         'GRN' = GRN_prob,
                         'CEACAM' = CEACAM_prob)
write.csv(list_of_datasets, "pathways_prob.csv")
