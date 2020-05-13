####################  插补性状缺失值(set_SNP_trait使用)【只需插补uncor_trait】  ###########
#########   输入值是性状名称（字符串），性状的pearson_R表格，作为解释变量的性状测量数据，以及插补时所需的表型测量数据

impute_pheNA <- function(trait,pearson_r,data,new_all_phenotype_measure){
  phe_row <- match(trait,rownames(pearson_r))

  high_cor <- ""
  value <- 0.9
  while (length(high_cor)<=1) {
    high_cor <- pearson_r[,which(pearson_r[phe_row,]>=value)] %>%
      colnames() %>%
      setdiff(.,trait)
    value <- value-0.01
  }

  high_cor_num  <- match(high_cor,colnames(new_all_phenotype_measure))
  high_cor_data <- new_all_phenotype_measure[,high_cor_num , with=FALSE]

  phe_data  <- new_all_phenotype_measure[, phe_row , with=FALSE]

  lm_data <- cbind(high_cor_data,phe_data)
  colnames(lm_data) <- c(colnames(high_cor_data),"phe")

  anova_imput <- rpart(phe~.,data=lm_data[!is.na(lm_data$phe),],na.action = na.omit,method = "anova")
  temp <- dim(lm_data)[2]-1
  anova_pre <- predict(anova_imput,lm_data[is.na(lm_data$phe),1:temp])

  phe_na <- which(is.na(phe_data))
  mean <- mean(unlist(phe_data),na.rm=TRUE)
  sd <- sd(unlist(phe_data),na.rm=TRUE)

  rep <- rnorm(3000,mean = 0,sd=sd)

  rep <- rep[rep>=min(phe_data,na.rm = TRUE)-mean& rep<= max(phe_data,na.rm = TRUE)-mean] %>%
    sample(.,length(phe_na),replace = FALSE)

  impute_phe <- anova_pre + rep

  data[is.na(data)] <- impute_phe
  return(data)
}


















