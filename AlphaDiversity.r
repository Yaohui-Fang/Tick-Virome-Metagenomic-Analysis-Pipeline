###α多样性计算
#清除所有变量
rm(list=ls())
#加载vegan包
library(vegan)
#读入物种数据
otu<-read.table('totalcp.txt',header = T,row.names = 1,check.names=F)
#Shannon 指数
Shannon<-diversity(otu, index = "shannon", MARGIN = 2, base = exp(1))
#Simpson 指数
Simpson<-diversity(otu, index = "simpson", MARGIN = 2, base = exp(1))
#Richness 指数
Richness <- specnumber(otu,MARGIN = 2)
#合并
index<-as.data.frame(cbind(Shannon,Simpson,Richness))
#接下来分析的多样性指数一般不作为重点分析对象，但既然要写，就整理的完整一些
#转置物种数据
totu<-t(otu)
totu<-ceiling(as.data.frame(t(otu)))
#多样性指数
obs_chao1_ace<-t(estimateR(totu))
obs_chao1_ace<-obs_chao1_ace[rownames(index),]
index$Chao1<-obs_chao1_ace[,2]
index$Ace<-obs_chao1_ace[,4]
index$Sobs<-obs_chao1_ace[,1]
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(otu == 1) / colSums(otu)
#合并、导出数据
write.table(cbind(sample=c(rownames(index)),index),'totalcp.alpha.txt',row.names = F,sep = '\t',quote = F)