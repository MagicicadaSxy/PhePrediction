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

######  ��ȡȫ����״�������ݣ�����ȱʧֵ�����������������
all_phenotype_measure <- fread("F:\\Magicicada\\thesis\\test\\mouse\\2016\\phenotypes\\CFW_measures.txt",header = T)
row_name <- unlist(all_phenotype_measure[,1])
attributes(row_name) <- NULL
rownames(all_phenotype_measure) <- row_name
all_phenotype_measure <- all_phenotype_measure[,-1] ####����С��ı�������

phe_cor <- corr.test(all_phenotype_measure,use = "pairwise",method = "pearson")  ###�������������
pearson_r <- phe_cor$r


###############################   ���Ŵ�������0.35����״ɸѡ������Ϊexemplary_trait
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

exemplar_trait <- all_phe[which(all_phe$Heritability>0.35),]  ### ���Ŵ�������0.35����״��Ϊ�о�����


############################  ѡ���ÿ��exemplary_trait��uncorrelated traits
trait_pos <- match(unlist(sort(exemplar_trait[,2])),colnames(all_phenotype_measure))

out <- list()
for (i in 1:27){
  m <- trait_pos[i]
  n <- pearson_r[m,-m]
  uncor_trait <- n[abs(n)<=0.30]     ### pearson_r�ľ���ֵС��0.3
  #uncor_trait <- names(uncor_trait)  ### �����Ҫ��Ӧ��r�����԰���һ��ע�͵�
  out[[i]] <- uncor_trait
}

names(out) <- unlist(sort(exemplar_trait[,2]))  ####��һ��27����б���ÿһ��Ϊһ��exemplay trait �� uncorrelated traits,����Ϊexemplary trait��

setwd("E:\\thesis")
save(out,pearson_r,file="traits_pearson_r.RData")