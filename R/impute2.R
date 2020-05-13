####################  插补性状缺失值(trait_to_trait使用)[exemplar和uncor_trait都需要插补]  ###########
#########   输入值是性状名称（字符串），性状的pearson_R表格，性状的全部测量数据

impute_pheNA2 <- function(trait,pearson_r,data){
  phe_row <- match(trait,colnames(all_phenotype_measure))
  
  high_cor <- ""
  value <- 0.9
  while (length(high_cor)<=2) {
    high_cor <- pearson_r[,which(pearson_r[phe_row-1,]>=value)] %>%
      colnames() %>%
      setdiff(.,trait)
    value <- value-0.01
  }
  
  high_cor_num  <- match(high_cor,colnames(all_phenotype_measure))
  high_cor_data <- all_phenotype_measure[,high_cor_num , with=FALSE]
  
  phe_data  <- all_phenotype_measure[, phe_row , with=FALSE]
  
  lm_data <- cbind(high_cor_data,phe_data)
  colnames(lm_data) <- c(colnames(high_cor_data),"phe")
  
  anova_imput <- rpart(phe~.,data=lm_data[!is.na(lm_data$phe),],na.action = na.omit,method = "anova")
  temp <- dim(lm_data)[2]-1
  anova_pre <- predict(anova_imput,lm_data[is.na(lm_data$phe),1:temp])
  
  phe_na <- which(is.na(phe_data))
  mean <- mean(unlist(phe_data),na.rm=TRUE)
  sd <- sd(unlist(phe_data),na.rm=TRUE)
  
  rep <- rnorm(30000,mean = 0,sd=sd)
  
  rep <- rep[rep>=min(phe_data,na.rm = TRUE)-mean& rep<= max(phe_data,na.rm = TRUE)-mean] %>%
    sample(.,length(phe_na),replace = FALSE)
  
  impute_phe <- anova_pre + rep
  
  data[is.na(data)] <- impute_phe
  return(data)
}








