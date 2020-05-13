###################  model_SNP+trait_cv  ##########

model_SNP_Trait_cv <- function(out_data){
  SNP_phe  <- round(as.data.frame(out_data),2)
  
  if(dim( SNP_phe )[2]==1){
    adjusted_R2 <- NA
    return(adjusted_R2)
  }
  
  col_num           <-  dim(SNP_phe)[2]-1
  SNP_data <- as.matrix(SNP_phe[,1:col_num])
  phe_data <- as.matrix(SNP_phe[,col_num+1])
  
  uni <- unique(phe_data)
  if(length(uni) == 1){
    adjusted_R2 <- NA
    return(adjusted_R2)
  }
  
  if(dim(SNP_data)[2] == 1){
    adjusted_R2 <- NA
    return(adjusted_R2)
  }
  
  if(FALSE){
    folds <- createFolds(SNP_phe$phenotype,k=10)
    adjusted_R2 <- lapply(folds, function(x){
      train <- as.matrix(SNP_phe[x,])
      test  <- as.matrix(SNP_phe[-x,])
      
      x_col <- dim(test)[2]-1
      train_x <- train[,1:x_col]
      train_y <- train[,x_col+1]
      test_x <- test[,1:x_col]
      test_y <- test[,x_col+1]
      
      model_cv <- cv.glmnet(train_x,train_y,family = "gaussian",type.measure = "mse",nfolds = 10,parallel = TRUE)
      
      train_R2 <- which(model_cv$lambda==model_cv$lambda.min) %>%
        model_cv$glmnet.fit$dev.ratio[.]
      
      #print(paste("train_R2:",train_R2))
      
      y_predict <- predict(model_cv,newx = test_x)
      sst <- sum((test_y-mean(test_y))^2)
      sse <- sum((y_predict-test_y)^2)
      test_R2 <- 1-sse/sst
      #print(paste("test_r2:",test_R2))
      
      out <- c(train_R2,test_R2)
      return(out)
    })
    
    train_R2 <- lapply(adjusted_R2, function(x){
      f0 <- x[1]
    }) %>%
      unlist() %>%
      max()
    
    test_R2 <- lapply(data, function(x){
      f0 <- x[2]
    }) %>%
      unlist() %>%
      max()
    
    print(paste("max_train_R2:",train_R2))
    print(paste("max_test_R2:",test_R2))
    
    new_R2 <- c(train_R2,test_R2)
    return(new_R2)
  } ### 十折交叉验证+test data验证 PS：test data带入模型后，效果普遍较差，输入为两元素的向量，第一个为train R2,第二个为test R2.
  
  if(TRUE){
    model_cv <- cv.glmnet(SNP_data,phe_data,family = "gaussian",type.measure = "mse",nfolds = 10,parallel = TRUE)
    
    adjusted_R2 <- which(model_cv$lambda==model_cv$lambda.min) %>%
      model_cv$glmnet.fit$dev.ratio[.] 
    
    return(adjusted_R2)
  } ### 仅十折交叉验证 ，输出为一个数字。
}





