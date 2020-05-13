library(dplyr)
model_traitTOtrait <- function(n){
  exemplar <- names(out)[n] %>%
    match(.,colnames(all_phenotype_measure)) %>%
    all_phenotype_measure[,.,with=FALSE] %>%
    as.matrix()  

  uncor_trait_data <- out[[n]] %>%
    match(.,colnames(all_phenotype_measure)) %>%
    all_phenotype_measure[,.,with=FALSE] %>%
    as.matrix()   


  trait    <-  names(out)[n]
  exemplar <- impute_pheNA2(trait = trait,pearson_r = pearson_r,data = exemplar)

  print(paste(trait,"h2;",exemplar_trait[n,13]))

  for(i in 1:dim(uncor_trait_data)[2]){
    phe <- colnames(uncor_trait_data)[i]
    uncor_trait_data[,i] <- impute_pheNA2(trait=phe,pearson_r = pearson_r,data = uncor_trait_data[,i])
  }

  if(TRUE){
    folds <- createFolds(exemplar,k=10)
    adjusted_R2 <- lapply(folds, function(x){
      train_y <- exemplar[-x,]
      test_y  <- exemplar[x,]
      train_x <- uncor_trait_data[-x,]
      test_x <- uncor_trait_data[x,]

      model_cv <- cv.glmnet(train_x,train_y,family = "gaussian",type.measure = "mse",nfolds = 10,parallel = TRUE)

      train_R2 <- which(model_cv$lambda==model_cv$lambda.min) %>%
        model_cv$glmnet.fit$dev.ratio[.]

      print(paste("train_R2:",train_R2))

      y_predict <- predict(model_cv,newx = test_x)
      sst <- sum((test_y-mean(test_y))^2)
      sse <- sum((y_predict-test_y)^2)
      test_R2 <- 1-sse/sst
      print(paste("test_r2:",test_R2))

      out <- c(train_R2,test_R2)
      return(out)
    })

  }  ?

  if(FALSE){
    model_cv <- cv.glmnet(uncor_trait_data,exemplar,type.measure = "mse",nfolds = 20,family = "gaussian",parallel = TRUE)

    y_predict <- predict(model_cv,newx = uncor_trait_data,s = "lambda.min")
    sst <- sum((exemplar-mean(exemplar))^2)
    sse <- sum((y_predict-exemplar)^2)
    R2 <- 1-sse/sst

    coef_num <- as.matrix(coef(model_cv,s="lambda.min"))
    coef_num <- coef_num[coef_num!=0]
    coef_num <- coef_num[-1]
    cv_adjusted_R2 <- 1- (sse/(dim(exemplar)[1]-length(coef_num)-1))/(sst/(dim(exemplar)[1]-1))
    print(paste("cv_r2:",cv_adjusted_R2))
  } 

  model_cv <- cv.glmnet(uncor_trait_data,exemplar,type.measure = "mse",nfolds = 20,family = "gaussian",parallel = TRUE)

  y_predict <- predict(model_cv,newx = uncor_trait_data,s = "lambda.min")
  sst <- sum((exemplar-mean(exemplar))^2)
  sse <- sum((y_predict-exemplar)^2)
  R2 <- 1-sse/sst

  coef_num <- as.matrix(coef(model_cv,s="lambda.min"))
  coef_num <- coef_num[coef_num!=0]
  coef_num <- coef_num[-1]
  cv_adjusted_R2 <- 1- (sse/(dim(exemplar)[1]-length(coef_num)-1))/(sst/(dim(exemplar)[1]-1))
  print(paste("cv_R2:",cv_adjusted_R2))

  lm_data <- as.data.frame(cbind(uncor_trait_data,exemplar))
  colnames(lm_data) <- c(colnames(uncor_trait_data),"phe")
  model_lm <- lm(phe~.,lm_data)
  lm_adjusted_R2 <- summary(model_lm)["adj.r.squared"] %>% unlist()
  attributes(lm_adjusted_R2) <- NULL
  print(paste("lm_R2:",lm_adjusted_R2))
  print("==================================================================================")

  a <- c(lm_adjusted_R2,cv_adjusted_R2)
  return(a)
}  ###??Ԫ???Իع??͵????????ع?ģ??
