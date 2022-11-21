rm(list = ls())
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/inferCNV")

##加载R包
library(infercnv)
library(AnnoProbe)

##导入数据
dat <- readRDS("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/seurat流程/Epi_exp.rds")
groupFiles <- read.table("groupFiles.txt", sep = '\t', row.names = 1)

##准备分析所需文件
geneInfor <- annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)

geneInfor <- geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]

length(unique(geneInfor[,1]))
head(geneInfor)

dat <- dat[rownames(dat) %in% geneInfor[,1],]
dat <- dat[match(geneInfor[,1], rownames(dat) ),] 
dim(dat)

expFile='expFile.txt'
write.table(dat,file = expFile,sep = "\t",quote = F)

head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = "\t",quote = F,col.names = F,row.names = F)

##创建inferCNV对象
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expFile,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = c("normal"))

##运行标准的inferCNV流程
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir = "./inferCNV_output",
                              cluster_by_groups = T,
                              denoise = T, #去噪
                              HMM = T) # 是否基于HMM预测CNV

