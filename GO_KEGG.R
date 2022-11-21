rm(list = ls()) 
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/GO_KEGG")

##加载R包
{
  library(AnnotationHub)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
}

##导入数据
geneSet1 <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/DNB/DNBgenelist.csv",
                     header = T)
geneSet2 <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/DNB/DNBneighbor.csv",
                     header = T)

##数据处理
geneSet1$ENSEMBL <- mapIds(org.Hs.eg.db,
                           geneSet1$gene,
                           "ENTREZID",
                           "SYMBOL")
geneSet1 <- na.omit(geneSet1)

geneSet2$ENSEMBL <- mapIds(org.Hs.eg.db,
                           geneSet2$gene,
                           "ENTREZID",
                           "SYMBOL")
geneSet2 <- na.omit(geneSet2)

##GO富集
ego<- enrichGO(geneSet1$ENSEMBL,
               OrgDb = 'org.Hs.eg.db',
               ont = 'BP',
               readable = T,
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH',
               minGSSize = 5,
               maxGSSize = 500,
               qvalueCutoff = 0.05)

go_fil <- ego@result
go_fil <- go_fil %>% subset(.,qvalue<0.05)
go_fil <- go_fil[order(go_fil$Count,decreasing = T),]
write.csv(go_fil,'DNB_BP_results.csv')

##GO结果可视化
go_fil2 <- read.csv("DNB_BP_select.csv", row.names = 1)
datbar <- go_fil2[,c(2,7,9)]
colnames(datbar) <- c("Term","FDR","Count")
datbar$Term <- factor(datbar$Term,levels = rev(datbar$Term))
plot1 <- ggplot(datbar,aes(x=Count,y=Term,fill=FDR))+
  geom_bar(stat = 'identity',mapping = aes(fill=FDR))+
  theme_test()+
  scale_colour_gradient(low="#d7dded",high="#5e79ba",aesthetics = "fill")+
  theme(axis.text = element_text(size = 10,angle = 0,colour = 'black'))+
  xlab('Count')+ylab('')
plot1      
ggsave('go_enrich.pdf',plot1,width = 10,height = 6)


##KEGG富集
# R.utils::setOption("clusterProfiler.download.method",'auto')

ekk<- enrichKEGG(geneSet2$ENSEMBL,
                 organism = "hsa",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'BH',
                 minGSSize = 5,
                 maxGSSize = 500,
                 qvalueCutoff = 0.05)

kegg_fil <- ekk@result
kegg_fil <- kegg_fil %>% subset(.,qvalue<0.05)
kegg_fil <- kegg_fil[order(kegg_fil$Count,decreasing = T),]
write.csv(kegg_fil,'DNBneighbor_KEGG_results.csv')

##KEGG结果可视化
kegg_fil2 <- read.csv("DNBneighbor_KEGG_select.csv", row.names = 1)
databar <- kegg_fil2[,c(2,7,9)]
colnames(databar) <- c("Term","FDR","Count")
databar$Term <- factor(databar$Term,levels = rev(databar$Term))
plot2 <- ggplot(databar,aes(x=Count,y=Term,fill=FDR))+
  geom_bar(stat = 'identity',mapping = aes(fill=FDR))+
  theme_test()+
  scale_colour_gradient(low="#d7dded",high="#5e79ba",aesthetics = "fill")+
  theme(axis.text = element_text(size = 10,angle = 0,colour = 'black'))+
  xlab('Count')+ylab('')
plot2      
ggsave('kegg_enrich.pdf',plot2,width = 10,height = 6)

