rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)
setwd("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/MFUZZ")

##加载R包
library(Mfuzz)

##导入数据
exp <- read.csv("D:/11/Subject/adenoma-carcinoma/GSE161277_Epi/monocle/EpiALL_DNBcount_raw(labelled_pseudotime).csv",
                header = T, row.names = 1)
DNB_neighbor <- read.csv("DNB_neighbor.csv")

index1 <- which(row.names(exp) %in% DNB_neighbor$neighbor)
exp_MFUZZ <- exp[index1,]
exp_MFUZZ <- rbind(exp[1,],exp_MFUZZ)

names(exp_MFUZZ) <- exp_MFUZZ[1,]
exp_MFUZZ <- exp_MFUZZ[-1,]
write.csv(exp_MFUZZ,"DNB_neighbor_countexp.csv")

##分群均值处理
index1 <- names(exp_MFUZZ) ==1
index2 <- names(exp_MFUZZ) ==2
index3 <- names(exp_MFUZZ) ==3
index4 <- names(exp_MFUZZ) ==4
index5 <- names(exp_MFUZZ) ==5
index6 <- names(exp_MFUZZ) ==6

exp_1 <- exp_MFUZZ[,index1]
exp_2 <- exp_MFUZZ[,index2]
exp_3 <- exp_MFUZZ[,index3]
exp_4 <- exp_MFUZZ[,index4]
exp_5 <- exp_MFUZZ[,index5]
exp_6 <- exp_MFUZZ[,index6]

Cluster1 <- apply(exp_1,1,mean)
Cluster2 <- apply(exp_2,1,mean)
Cluster3 <- apply(exp_3,1,mean)
Cluster4 <- apply(exp_4,1,mean)
Cluster5 <- apply(exp_5,1,mean)
Cluster6 <- apply(exp_6,1,mean)

MFUZZ_mean_exp <- data.frame(Cluster1,Cluster2,Cluster3,Cluster4,
                             Cluster5,Cluster6)
write.csv(MFUZZ_mean_exp,"DNB_neighbor_MFUZZ_mean_exp.csv")

##Mfuzz软聚类 
###构建对象
MFUZZ_mean_exp <- as.matrix(MFUZZ_mean_exp)
dat <- ExpressionSet(assayData = MFUZZ_mean_exp)

###排除超过25%的测量缺失的基因
dat.r <- filter.NA(dat, thres=0.25)

###用相应基因的平均值表达值替换剩余的缺失值
dat.f <- fill.NA(dat,mode="mean")
tmp <- filter.std(dat.f,min.std=0)

###标准化
dat.s <- standardise(tmp)

###聚类
set.seed(2022)

###手动定义聚类个数c
cl <- mfuzz(dat.s,c=16,m=1.25)

###作图
mfuzz.plot(dat.s,cl=cl,mfrow=c(4,4))

##保存MFUZZ对象
saveRDS(cl,"cl_v1.rds")
cl <- readRDS("cl_v1.rds")

##保存标准化的文件
a <- dat.s@assayData$exprs
write.csv(a,"exp_normalized_MFUZZ.csv")
