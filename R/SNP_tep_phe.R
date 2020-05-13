###### SNP、tep和phe数据提取，参数为m（性状序号 4:205）和SNP_class（high_SNP,random_SNP,random_SNP2）
###### 输出结果供SNP_dos_base函数使用

SNP_tep_phe <- function(m,SNPclass,setwd){
  setwd(setwd)

  pheName <- dir()                     #### 某表型与SNP相关性 前4个不是所需文件

  phe0 <- str_sub(pheName[m],1,-13)  #####行坐标表示某性状
  print(paste("load trait",phe0))

  h2 <- h2_data[which(h2_data[,3]==phe0),14]
  print(paste("h2:",h2))

  pheName <- paste(phe0,".mapping.txt",sep = "")

  a       <- fread(pheName,header = F,quote = "\"",col.names = c("V","chr","bp","-logP"))[,2:4]  ##读取某一性状与SNP的P值.
  phe     <- select(all_phenotype_measure,"Animal_ID",phe0)
  phe     <- na.omit(phe)  ####去除NA的小鼠

  print(paste("trait",phe0,"has",dim(phe)[1],"mice",sep=""))

  SNP6       <- a[which(a[,3]>6)]
  SNP_number <- dim(SNP6)[1]

  if(SNP_number==0){
    SNP_number <- round(dim(phe)[1]/3)
  }

  if(SNPclass == "high_SNP"){
    SNP <- SNP6
  }else if(SNPclass == "random_SNP"){
    SNP <- a[sample(dim(a)[1],SNP_number),]  #### 随机抽取的SNP数量等于高相关性的SNP数量
  }else if(SNPclass == "random_SNP2"){
    SNP <- a[sample(dim(a)[1],2*SNP_number),]
  }

  print(paste("total",dim(SNP)[1],"SNP",sep=""))

  rm(a)
  rm(SNP6)
  tep <- as.data.frame(table(SNP[,1]))  ####  检查SNP集中在哪些染色体[PS:染色体号是字符串]
  tep <- tep[,1]

  a <- str_c(tep,collapse = ",")
  print(paste("SNP on chr",a,sep=""))

  out <- list(SNP,tep,phe)

}



