model_traitTOtrait_svr <- function(n){
  exemplar <- names(out)[n] %>%
    match(.,colnames(all_phenotype_measure)) %>%
    all_phenotype_measure[,.,with=FALSE] %>%
    as.matrix()  

  uncor_trait_data <- out[[n]] %>%
    match(.,colnames(all_phenotype_measure)) %>%
    all_phenotype_measure[,.,with=FALSE] %>%
    as.matrix()  


  trait <-  names(out)[n]
  exemplar <- impute_pheNA2(trait = trait,pearson_r = pearson_r,data = exemplar)

  for(i in 1:dim(uncor_trait_data)[2]){
    phe <- colnames(uncor_trait_data)[i]
    uncor_trait_data[,i] <- impute_pheNA2(trait=phe,pearson_r = pearson_r,data = uncor_trait_data[,i])
  }
  SNP_phe <- as.data.frame(cbind(uncor_trait_data,exemplar))
  col_num           <-  dim(SNP_phe)[2]-1
  SNP_name <- paste("trait",c(1:col_num),sep = "")
  colnames(SNP_phe) <-  c(SNP_name,"phenotype")
  SNP_phe <- as.data.frame(SNP_phe)

  test_num <- dim(SNP_phe)[1]
  test_mice <- sample(1:test_num,test_num/5)

  train <- SNP_phe[-test_mice,]
  test_x  <- SNP_phe[test_mice,1:col_num]
  test_y  <- SNP_phe[test_mice,col_num+1]
  model_svm <- ksvm(phenotype~.,train,kernel= "rbfdot",type="eps-bsvr")

  train_R2 <- 1-model_svm@error
  print(paste("train_R2:",train_R2))

  y_predict <- predict(model_svm, test_x)
  sst <- sum((test_y-mean(test_y))^2)
  sse <- sum((y_predict-test_y)^2)
  test_R2 <- 1-sse/sst
  print(paste("test_r2:",test_R2))

  R2 <- c(train_R2,test_R2)
  return(R2)
}