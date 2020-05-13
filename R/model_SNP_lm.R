#####  利用SNP剂量数据（SNP_dos_base函数输出）构建lm模型
#### 拟合生成的多元线性模型无法用于test data。

model_SNP_lm <- function(out_data){
  myData2 <- as.data.frame(out_data[2])

  if(dim(myData2)[2]==1){
    adjusted_R2 <- NA
    return(adjusted_R2)
  }

  SNP_phe           <-  t(round(myData2,1))

  col_num           <-  dim(SNP_phe)[2]-1
  col_name          <-  paste("SNP",c(1:col_num),sep = "")
  colnames(SNP_phe) <-  c(col_name,"phetype")

  SNP_phe <- as.data.frame(SNP_phe)
  model_lm <- lm(phetype ~.,SNP_phe)

  adjusted_R2 <- summary(model_lm)["adj.r.squared"] %>% unlist()
  attributes(adjusted_R2) <- NULL

  return(adjusted_R2)
}



