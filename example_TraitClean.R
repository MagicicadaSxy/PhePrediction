library(psych)
library(openxlsx)
library(data.table)
library(dplyr)
library(stringr)
library(glmnet)
library(rpart)

library(doParallel)
library(kernlab)
library(caret)

######  读取全部性状测量数据，处理缺失值，并两两计算相关性
all_phenotype_measure <- fread("F:\\Magicicada\\thesis\\test\\mouse\\2016\\phenotypes\\CFW_measures.txt",header = T)
row_name <- unlist(all_phenotype_measure[,1])
attributes(row_name) <- NULL
rownames(all_phenotype_measure) <- row_name
all_phenotype_measure <- all_phenotype_measure[,-1] ####所有小鼠的表型数据

phe_cor <- corr.test(all_phenotype_measure,use = "pairwise",method = "pearson")  ###两两计算相关性
pearson_r <- phe_cor$r


###############################   将遗传力大于0.35的性状筛选出来作为exemplary_trait
all_phe <- read.xlsx("E:\\R-project\\test\\H2.xlsx",)
all_phe <- all_phe[3:203,c(1,3:14)]
colnames(all_phe) <- all_phe[1,]
all_phe <- all_phe[-1,]
all_phe[,c(3:10,13)] <-lapply(all_phe[,c(3:10,13)],as.numeric)

for (i in 2:200) {
  if (is.na(all_phe[i,1])){
    all_phe[i,1] <- all_phe[i-1,1]
  }
}

exemplar_trait <- all_phe[which(all_phe$Heritability>0.35),]  ### 将遗传力大于0.35的性状作为研究对象


############################  选择出每个exemplary_trait的uncorrelated traits
trait_pos <- match(unlist(sort(exemplar_trait[,2])),colnames(all_phenotype_measure))

out <- list()
for (i in 1:27){
  m <- trait_pos[i]
  n <- pearson_r[m,-m]
  uncor_trait <- n[abs(n)<=0.30]     ### pearson_r的绝对值小于0.3
  #uncor_trait <- names(uncor_trait)  ### 如果需要对应的r，可以把这一行注释掉
  out[[i]] <- uncor_trait
}

names(out) <- unlist(sort(exemplar_trait[,2]))  ####是一个27层的列表，每一层为一个exemplay trait 的 uncorrelated traits,层名为exemplary trait。

setwd("E:\\thesis")
save(out,pearson_r,file="traits_pearson_r.RData")
