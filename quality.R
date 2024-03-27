setwd("D://LILab//SC蛋白组//231122//final_result")

data <- read.table("ev_updated.txt", header = T,
                   stringsAsFactors = F, sep = "\t")
t = colnames(data)
RI = paste0('RI',seq(1,16))
t[grep('Reporter.intensity.corrected',colnames(data))] = RI
colnames(data) = t
############# 2. 加载所需R包
suppressMessages(library(scp))
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(impute))
suppressMessages(library(sva))
suppressMessages(library(scater))
suppressMessages(library(SCP.replication))
############# 3. 注释16个chanel对应细胞类型
SampleAnnotation = readRDS('data/sampleannotation.rds')
############# 4. 构建scp project
scp <- readSCP(featureData = data,
               colData = SampleAnnotation,
               channelCol = "Channel",
               batchCol = "Raw.file",
               removeEmptyCols = TRUE)