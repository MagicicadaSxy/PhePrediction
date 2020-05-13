#####  利用SNP剂量数据及其uncor_trait数据（set_SNP_trait函数输出）构建lm模型
#### 拟合生成的多元线性模型无法用于test data。

model_SNP_Trait_lm <- function(out_data){
  SNP_phe  <- round(as.data.frame(out_data),2)
  
  if(dim( SNP_phe )[2]==1){
    adjusted_R2 <- NA
    return(adjusted_R2)
  }
  
  col_num                    <-  dim(SNP_phe)[2]
  colnames(SNP_phe)[col_num] <-  c("phenotype")
  
  SNP_phe  <- as.data.frame(SNP_phe)
  model_lm <- lm(phenotype ~.,SNP_phe)
  
  adjusted_R2 <- summary(model_lm)["adj.r.squared"] %>% unlist()
  attributes(adjusted_R2) <- NULL
  
  return(adjusted_R2)
}
