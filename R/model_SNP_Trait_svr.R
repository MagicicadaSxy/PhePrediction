#####  利用SNP剂量数据及其uncor_trait数据（set_SNP_trait函数输出）构建svr模型


model_SNP_Trait_svr <- function(out_data){
  SNP_phe  <- round(as.data.frame(out_data),2)

  if(dim( SNP_phe )[2]==1){
    adjusted_R2 <- NA
    return(adjusted_R2)
  }

  col_num           <-  dim(SNP_phe)[2]-1
  colnames(SNP_phe)[col_num+1] <-  c("phenotype")

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







