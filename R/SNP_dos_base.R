###### 传入项为SNP_tep_phe函数的输出结果（包含tep，SNP和phe），输出为一个列表，列表共两层，第一层为myData1，第二层为myData2

SNP_dos_base <- function(input){
  SNP <- as.data.frame(input[1])
  tep <- unlist(input[2])
  phe <- as.data.frame(input[3])

  if (length(tep)==0){
    out_data <- list(myData1 = 0,myData2 = 0)
    return(out_data)
  }

  setwd("F:\\Magicicada\\thesis\\test\\mouse\\2016\\Genotype Dosages\\dosages")
  SNP_base <- data.frame() ###每一列为一个小鼠样本，每一行为一个SNP样本
  SNP_dos  <- data.frame()

  for (i in tep) {
    temp <- paste("chr",i,".RData",sep = "")

    print(paste("chr",i,"loading"))

    load(temp)
    pruned_dosages <- round(pruned_dosages,2)
    SNP_pos   <- SNP[which(unlist(SNP[,1])==i),2]      
    row_num   <- na.omit(match(unlist(SNP_pos),pruned_pos[,2]))  ##匹配得到的行名 行名代表SNP
    SNP_base0 <- pruned_pos[row_num,] ##SNP的碱基
    SNP_base  <- rbind(SNP_base,SNP_base0)   ###SNP_base文件即为需要的基因型文件

    all_mice <- dos_name(nameList)
    emm      <- unlist(as.vector(phe[,1]))
    mice_num <- na.omit(match(emm,all_mice))
    SNP_dos0 <- pruned_dosages[row_num,mice_num]  ### 有表型数据的小鼠及其所需SNP的剂量

    if(is.null(dim(SNP_dos0))==TRUE){
      SNP_dos0 <- matrix(SNP_dos0,nrow = 1)
    }else{
      SNP_dos0 <- matrix(SNP_dos0,nrow = dim(SNP_dos0)[1])
    }

    SNP_dos <- rbind(SNP_dos,SNP_dos0)

    print(paste(i,"chr over"))
  }

  rm(SNP_base0)
  rm(SNP_dos0)
  rm(pruned_dosages)

  noPhe <- setdiff(emm,all_mice)
  noPhe <- match(noPhe,unlist(as.vector(phe[,1])))

  if(length(noPhe)==0){
    final_phe <- t(phe)
  }else{
    final_phe <- t(phe[-noPhe,])
  }

  colnames(SNP_dos) <-  final_phe[1,]

  myData1 <- cbind(SNP_base,SNP_dos) ###每一列为一个小鼠样本，每一行为一个SNP样本,无表型,前四列为染色体号，碱基位置，REF和ALT
  myData2 <- rbind(SNP_dos,as.numeric(final_phe[2,]))   #### 每一列为一个小鼠样本，每一行为一个SNP样本,无碱基,最后一行表型

  out_data <- list(myData1,myData2)
}






