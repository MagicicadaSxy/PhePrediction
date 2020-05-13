
###将剂量文件中的小鼠代号与表型文件中的小鼠代号统一

dos_name <- function(nameList){
  name_chr <- unlist(nameList)
  name_chr <- strsplit(name_chr,"_recal")
  name_chr <- sapply(name_chr,"[",1)
  name_chr <- gsub("___","/",name_chr)
}
