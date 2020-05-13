#################  SNP_trait  ####################
###### 传入项为m 和 SNP_dos_base函数的输出结果，输出为一个数据框,每行为一个mice，每一列：前半部分为uncor_trait，后半部分为SNP，最后一列为预测表型

set_SNP_trait <- function(m,out_data){
  new_all_phenotype_measure <- fread("F:\\Magicicada\\thesis\\test\\mouse\\2016\\phenotypes\\CFW_measures.txt",header = T)
  
  SNP_data <- as.data.frame(out_data[2])

  exemplar_name <- colnames(SNP_data)

  new_all_phenotype_measure <-match(exemplar_name,unlist(new_all_phenotype_measure[,1])) %>%
    new_all_phenotype_measure[.,]
  rownames(new_all_phenotype_measure) <- unlist(new_all_phenotype_measure[,1])
  new_all_phenotype_measure <- new_all_phenotype_measure[,-1]


  uncor_trait_data <- out[[m]] %>%
    match(.,colnames(new_all_phenotype_measure)) %>%
    new_all_phenotype_measure[,.,with=FALSE] %>%
    as.matrix()   #### 自变量数据
  
  ### 填补空缺值
  for(i in 1:dim(uncor_trait_data)[2]){
    phe <- colnames(uncor_trait_data)[i]
    uncor_trait_data[,i] <- impute_pheNA(trait=phe,pearson_r = pearson_r,data = uncor_trait_data[,i],new_all_phenotype_measure=new_all_phenotype_measure)
  }
  print(dim(uncor_trait_data))


  out <- cbind(uncor_trait_data,t(SNP_data))

  return(out)
}
